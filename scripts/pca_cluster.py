#!/usr/bin/env python
# From stdlib
import sys
import os
import argparse
from collections import Counter
from multiprocessing import Pool

# Generic third party libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import optuna
from optuna.trial import TrialState
import joblib

# Dimension reduction and clustering libraries
from umap import umap_ as umap
import hdbscan
import sklearn.cluster as cluster
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score, silhouette_score
from sklearn.feature_selection import mutual_info_classif
from sklearn.metrics import adjusted_mutual_info_score, normalized_mutual_info_score
from sklearn.feature_selection import VarianceThreshold

# Trajectory reading and analysis
import MDAnalysis
from MDAnalysis.analysis import rms


def make_parser():
    parser = argparse.ArgumentParser(description="Script to cluster protein structures using PCA and HDBSCAN")
    parser.add_argument("-g", "--gro", default=None, help="GRO file")
    parser.add_argument("-t", "--trr", default=None, help="TRR or XTC file")
    parser.add_argument("-d", "--dataset", default=None, help="dataset name")
    parser.add_argument("-n", "--n_trials", default=50, type=int, help="Number of optimization trials")
    return parser


def main(argv):
    parser = make_parser()
    args = parser.parse_args(argv[1:]) 

    basename = args.trr.replace(".trr", "").replace(".xtc", "").replace(".pdb", "")

    u = MDAnalysis.Universe(args.gro, args.trr, in_memory=True)  # always start with a Universe

    data = u.trajectory.coordinate_array 
    data = data.reshape(data.shape[0], data.shape[1] * 3)

    cluster(data, u, args.dataset, args.n_trials)


def calc_rmsd(u, cluster):

    u2 = u.copy()

    # select CA
    ca = u.select_atoms('name CA')
    ca2 = u2.select_atoms('name CA')
    n_ca = len(ca)

    #import math
    #nn = len(cluster)
    #print("perm", math.factorial(nn) / (math.factorial(2) * math.factorial(nn - 2)))
    data = []
    traj1 = u.trajectory[cluster]
    traj2 = u2.trajectory[cluster]
    for ts in traj1:     # iterate through all frames in u1
        for ts2 in traj2:     # iterate through all frames in u2
            if ts.frame < ts2.frame:
                #print(ts.frame, ts2.frame)
                r = rms.rmsd(ca.positions,  # coordinates to align
                    ca2.positions,  # reference coordinates
                    center=True,  # subtract the center of geometry
                    superposition=True)  # superimpose coordinates
                data.append(r)  # save RMSD value
    return np.mean(data)


def variance_threshold_selector(data, threshold=0.5):
    selector = VarianceThreshold(threshold)
    selector.fit(data)
    return data[data.columns[selector.get_support(indices=True)]]


def run_experiment(data, universe, n_components, min_cluster_size):
    clusterable_embedding = PCA(n_components=n_components).fit_transform(data)
    clusterer = hdbscan.HDBSCAN(
        min_samples=5,
        min_cluster_size=min_cluster_size,
        gen_min_span_tree=True,
    )
    clusterer.fit(clusterable_embedding)
    labels = clusterer.labels_
    clustered = (labels >= 0)
    rmsd = []
    for cluster in np.sort(np.unique(labels)):
        rmsd.append(float(calc_rmsd(universe, (labels == cluster))))

    mean_rmsd = np.mean(rmsd)
    d = {"relative_validity": float(clusterer.relative_validity_),
         "rmsd": list(rmsd),
         "mean_rmsd": float(mean_rmsd),
         "labels": [int(i) for i in labels],
         #"clusterable_embedding": clusterable_embedding,
         "cluster_persistence": list(clusterer.cluster_persistence_),
         "number_of_clusters": int(np.unique(labels).shape[0]),
         "frac_clustered": float(np.sum(clustered) / data.shape[0]),
         "silhouette_value": float(silhouette_score(clusterable_embedding[clustered, :], labels[clustered])),
         }
    #print(d)
    return d


class Objective:

    def __init__(self, data, universe):
        self.data = data
        self.universe = universe

    def __call__(self, trial):
        n_components = trial.suggest_categorical(
            "n_components", list(range(2, 51, 1)))
        min_cluster_size = trial.suggest_categorical(
            "min_cluster_size", list(range(100, 1100, 100)))
        d = run_experiment(self.data, self.universe, n_components,
                           min_cluster_size)
        for metric in d:
            print(metric, d[metric], type(metric), type(d[metric]))
            trial.set_user_attr(metric, d[metric])
        return d["mean_rmsd"]



def cluster(data, universe, dataset, n_trials=50):
    print(dataset)
    study_name = dataset.replace(".csv.gz", "") + "_pca.study"
    study = optuna.create_study(direction="minimize", study_name=study_name,
                                storage=f'sqlite:///{dataset}_pca.db', load_if_exists=True)
    if study.user_attrs.get("dataset", None) is None:
        study.set_user_attr("dataset", dataset)
    #data = pd.read_csv(dataset, index_col=0)
    objective = Objective(data, universe)
    try:
        n_complete_trials = len(
            [trial for trial in study.get_trials() if trial.state == TrialState.COMPLETE])
    except:
        n_complete_trials = 0
    print(f"Study {study_name} has {n_complete_trials} complete trials.")
    if n_complete_trials < n_trials:
        study.optimize(objective, n_trials=n_trials)
    joblib.dump(study, study_name)
    df = study.trials_dataframe()
    print(df.head())
    df.to_csv(f"{dataset}_pca_clustering.csv", index=False)

    df = df.sort_values(by="value", ascending=False)

    n_components = int(df.iloc[1, :]["params_n_components"])
    min_cluster_size = int(df.iloc[0, :]["params_min_cluster_size"])
    d = run_experiment(data, universe, n_components, min_cluster_size)
    output = pd.DataFrame(d["clusterable_embedding"], index=data.index)
    output["labels"] = d["labels"]
    output.to_csv(f"{dataset}_pca_embeddings.csv")
    print(d["frac_clustered"], d["silhouette_value"], d["relative_validity"])
    

if __name__ == '__main__':
    import sys
    main(sys.argv)
