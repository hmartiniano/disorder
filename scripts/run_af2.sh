#!/bin/bash
for i in $*
do
  colabfold_batch $i initial_structures/$(basename $i .fasta) --cpu
done

