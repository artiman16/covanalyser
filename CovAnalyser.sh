#!/bin/bash

cat Samples/*.fasta > samples.fasta
echo "Wait for alignment..."
mafft --quiet --reorder --6merpair --keeplength --addfragments samples.fasta reference.fa  > alignment.fa 
echo "Alignment completed!"
echo "Wait, analysis is running..."
python3 analysis.py
mv *.csv Results/
mv "alignment.fa" Results/
echo "Analysis comleted!"
