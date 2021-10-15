#!/bin/bash
for pdb in initial_structures/*pdb
do
  name=$(basename $pdb .pdb)
  mkdir -p $name
  cp -r Snakefile mdp/* $name 
  cp $pdb $name/protein.pdb
  cd $name 
  snakemake -s Snakefile --cores 8 --use-singularity >& ${name}.log
  cd ..
done
