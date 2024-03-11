#!/bin/bash

cat Samples/*.fasta > samples.fasta
echo "Подождите, идёт выравнивание..."
mafft --quiet --reorder --6merpair --keeplength --addfragments samples.fasta reference.fa  > alignment.fa 
echo "Выравнивание завершено!"
echo "Подождите, идёт анализ..."
python3 analysis.py
mv *.csv Results/
mv "alignment.fa" Results/
echo "Анализ завершён!"
