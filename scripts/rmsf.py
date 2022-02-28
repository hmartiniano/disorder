#!/usr/bin/env python
import os
import argparse

import numpy.linalg
import pandas as pd
import numpy as np
import networkx as nx
from networkx.algorithms.core import core_number

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis import distances

from MDAnalysis.analysis import align
import MDAnalysis.transformations as trans
from MDAnalysis.analysis.rms import RMSF

CUTOFF = 4.5  # CA-CA distance cutoff for defining residue contacts
#CUTOFF = 6  # CA-CA distance cutoff for defining residue contacts


def make_parser():
    parser = argparse.ArgumentParser(description="Script to create a Protein Strcture Network from inter-residue contacts")
    parser.add_argument("-g", "--gro", default=None, help="GRO file")
    parser.add_argument("-t", "--trr", default=None, help="TRR file")
    return parser


def main(argv):
    parser = make_parser()
    args = parser.parse_args(argv[1:]) 

    basename = args.trr.replace(".trr", "").replace(".xtc", "").replace(".pdb", "")

    u = mda.Universe(args.gro, args.trr)  # always start with a Universe

    NAC = 'backbone and resid 61-95'
    NTERM = 'backbone and resid 1-60'
    CTERM = 'backbone and resid 96-140'
	
    bb = u.select_atoms('protein and backbone')  # a selection (AtomGroup)
    nac = u.select_atoms(NAC) 
    nterm = u.select_atoms(NTERM)
    cterm = u.select_atoms(CTERM) 

    # select CA
    ca = u.select_atoms('name CA')
    n_ca = len(ca)

    protein = u.select_atoms("protein")

    # 1) the current trajectory contains a protein split across
    #    periodic boundaries, so we first make the protein whole and
    #    center it in the box using on-the-fly transformations

    #not_protein = u.select_atoms('not protein')
    #transforms = [trans.unwrap(protein),
    #              trans.center_in_box(protein, wrap=True),
    #              trans.wrap(not_protein)]
    #u.trajectory.add_transformations(*transforms)

    # 2) fit to the initial frame to get a better average structure
    #    (the trajectory is changed in memory)
    prealigner = align.AlignTraj(u, u, select="protein and name CA",
                                 in_memory=True).run()

    # 3) reference = average structure
    ref_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
    # make a reference structure (need to reshape into a 1-frame
    # "trajectory")
    reference = mda.Merge(protein).load_new(ref_coordinates[:, None, :],
                                        order="afc")

    aligner = align.AlignTraj(u, reference,
                          select="protein and name CA",
                          in_memory=True).run()

    calphas = protein.select_atoms("name CA")
    rmsfer = RMSF(calphas, verbose=True).run()


    data = pd.DataFrame(rmsfer.results.rmsf, index=calphas.resnums)
    data.to_csv("rmsfer.csv")

if __name__ == '__main__':
    import sys
    main(sys.argv)
