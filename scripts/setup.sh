#!/bin/bash
Help()
{
   # Display Help
   echo "Script to setup a series of molecular dynamics simulations from PDB files."
   echo
   echo "Syntax: run.sh [-i|s|w|h|v|V]"
   echo "options:"
   echo "i     directory with AF2 output directories."
   echo "s     Snakefile."
   echo "w     workdir."
   echo "h     Print this Help."
   echo "v     Verbose mode."
   echo "V     Print software version and exit."
   echo
}

while getopts p:s:w:h:v:V flag
do
    case "${flag}" in
        i) input_dir=${OPTARG};;
        s) snakefile=${OPTARG};;
        w) workdir=${OPTARG};;
        h) # Display help and exit
           Help
           exit;;
       \?) # Invalid option
           echo "Error: Invalid option"
           exit;;
    esac
done

input_dir=${workdir:-initial_structures}
snakefile=${snakefile:-workflows/md/Snakefile}
workdir=${workdir:-workdir}
echo $workdir
mkdir -p $workdir

for d in $(ls $input_dir)
do
  mkdir -p ${workdir}/${d}
  cp -r ${snakefile} ${workdir}/${d}
  cp -r ${input_dir}/${d}/${d}_unrelaxed_rank_1_model_1.pdb ${workdir}/${d}/protein.pdb
  cp -r $(dirname ${snakefile})/mdp/*mdp ${workdir}/${d}
done
