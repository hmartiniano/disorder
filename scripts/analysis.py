#!/usr/bin/env python
import argparse

import numpy.linalg
import pandas as pd

import MDAnalysis
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import contacts



def make_parser():
    parser = argparse.ArgumentParser(description="Script to calculate Rgyr, RMSD, Native Contacts and End-to-End distance")
    parser.add_argument("-g", "--gro", default=None, help="GRO file")
    parser.add_argument("-t", "--trr", default=None, help="TRR file")
    return parser


def main(argv):
    parser = make_parser()
    args = parser.parse_args(argv[1:]) 

    u = MDAnalysis.Universe(args.gro, args.trr)  # always start with a Universe

    NAC = 'backbone and resid 61-95'
    NTERM = 'backbone and resid 1-60'
    CTERM = 'backbone and resid 96-140'
	
    # RMSD
    R = rms.RMSD(u,  # universe to align
		 u,  # reference universe or atomgroup
		 select='backbone',  # group to superimpose and calculate RMSD
		 groupselections=[NAC, NTERM, CTERM],  # groups for RMSD
		 ref_frame=0)  # frame index of the reference
    R.run()
    rmsd_df = pd.DataFrame(R.rmsd,
                  columns=['Frame', 'Time (ns)',
                           'Backbone', 'NAC',
                           'NTERM', 'CTERM'])
    rmsd_df.to_csv("rmsd.csv", index=False)

    # Native Contacts
    q1q2 = contacts.q1q2(u, 'name CA', radius=6).run()

    q1q2_df = pd.DataFrame(q1q2.timeseries,
                       columns=['Frame',
                                'Q1',
                                'Q2'])
    q1q2_df.to_csv("contacts.csv", index=False)

    # can access via segid (4AKE) and atom name
    # we take the first atom named N and the last atom named C
    nterm = u.select_atoms('protein and name N')[0]
    cterm = u.select_atoms('protein and name C')[-1]

    bb = u.select_atoms('protein and backbone')  # a selection (AtomGroup)
    data = []
    for ts in u.trajectory:     # iterate through all frames
        r = cterm.position - nterm.position # end-to-end vector from atom positions
        d = numpy.linalg.norm(r)  # end-to-end distance
        rgyr = bb.radius_of_gyration()  # method of AtomGroup
        print("frame = {0}: d = {1} A, Rgyr = {2} A".format(
              ts.frame, d, rgyr))
        data.append((ts.frame, d, rgyr, 
            NAC.radius_of_gyration(),  
            NTERM.radius_of_gyration(),  
            CTERM.radius_of_gyration(),  
            ))
    data = pd.DataFrame(data, names=["Frame", "e-e distance", "Backbone Rgyr", "NAC Rgyr", "NTERM Rgyr", "CTERM Rgyr"])
    data = pd.concat((rmsd_df, q1q2_df.drop(columns="Frame"), data.drop(columns="Frame")), axis=1)
    data.to_csv("analysis.csv", index=False)
if __name__ == '__main__':
    import sys
    main(sys.argv)
