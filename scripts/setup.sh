#!/bin/bash
Help()
{
   # Display Help
   echo "Script to setup a series of molecular dynamics simulations from PDB files."
   echo
   echo "Syntax: run.sh [-p|s|w|h|v|V]"
   echo "options:"
   echo "p     PDB files."
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
        p) pdb_files=${OPTARG};;
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

workdir=${workdir:-wordkdir}
echo $workdir
mkdir -p $workdir

for f in $(ls $pdb_files/*.pdb)
do
  mkdir -p ${workdir}/${i}
  cp -r ${snakefile} ${workdir}/${i}
  cp -r ${f} ${workdir}/${i}/protein.pdb
  cp -r $(dirname ${snakefile})/mdp ${workdir}/${i}
done
