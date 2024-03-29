Disorder
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/hmartiniano/disorder/workflows/CI/badge.svg)](https://github.com/hmartiniano/disorder/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/hmartiniano/disorder/branch/master/graph/badge.svg)](https://codecov.io/gh/hmartiniano/disorder/branch/master)


A set of tools and pipelines for molecular dynamics simulations of protein systems.
This is mainly aimed at simulating instrinsically disordered proteins, though it is general enough to be used for any protein class.


Installation
============

First install conda, if it is not installed already.

Clone the repo from github:

```
git clone https://github.com/hmartiniano/Disorder
```

Create conda environment

```
cd disorder
conda env create -f environment.yml
```

Usage
=====

0. Activate the conda environment:
```
conda activate -n disorder 
```

1. Edit `notebooks/explore_variants.ipynb` and input a gene of interest and, optionally, known mutations.

2. Run the structure prediction pipeline:
```
./scripts/run_af2.sh <fasta_files> 
```

3. Run the system setup for a given pipeline (atomistic or coarse-grained):
```
snakemake -s workflows/setup/md/Snakefile.aa 
```

4. cd to the working directory and use snakemake to run individual simulations:
```
cd $workdir 
cd $mutation
snakemake -c <number of cores>
```

5. Run the analysis scripts in the ```scripts``` directory:


### Copyright

Copyright (c) 2022, Hugo Martiniano, Nuno Galamba


#### Acknowledgements
 
This software was produced wth support from [BioISI](htts://bioisi.pt).
