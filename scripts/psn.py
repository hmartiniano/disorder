#!/usr/bin/env python
import os
import argparse

import numpy.linalg
import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms.core import core_number

import MDAnalysis
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis import distances

CUTOFF = 4.5  # CA-CA distance cutoff for defining residue contacts
#CUTOFF = 6  # CA-CA distance cutoff for defining residue contacts


def make_parser():
    parser = argparse.ArgumentParser(description="Script to create a Protein Structure Network from inter-residue contacts")
    parser.add_argument("-g", "--gro", default=None, help="GRO file")
    parser.add_argument("-t", "--trr", default=None, help="TRR file")
    return parser


def main(argv):
    parser = make_parser()
    args = parser.parse_args(argv[1:]) 

    basename = args.trr.replace(".trr", "").replace(".xtc", "").replace(".pdb", "")

    u = MDAnalysis.Universe(args.gro, args.trr)  # always start with a Universe

    # select CA
    ca = u.select_atoms('name CA')
    n_ca = len(ca)

    data = []
    output_path='contact_matrix.csv'
    for n, ts in enumerate(u.trajectory):     # iterate through all frames
        dist = contacts.distance_array(ca.positions, ca.positions)
        cm = contacts.contact_matrix(dist, CUTOFF)
        pd.DataFrame(cm.astype(int)).to_csv(output_path, mode='a', header=not os.path.exists(output_path))
        psn = nx.from_numpy_array(cm)
        psn.remove_edges_from(nx.selfloop_edges(psn))
        #print(psn.number_of_nodes(), psn.number_of_edges())
        c = core_number(psn)
        data.append([c[i] for i in range(n_ca)])
        assert len(c) == n_ca
    data = pd.DataFrame(data)
    data.to_csv("core_number.csv")

if __name__ == '__main__':
    import sys
    main(sys.argv)
