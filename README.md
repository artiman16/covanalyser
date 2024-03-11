# CovAnalyser: SARS-CoV-2 genotyping program

## Introduction

CovAnalyser is a software tool to find mutations in SARS-CoV-2 Spike genes and determine genetic lineages of isolates.

## Installation (Linux)
```
git clone https://github.com/artiman16/covanalyser.git
conda env create -f environment.yml
conda activate covanalyser
<<<<<<< HEAD
=======

>>>>>>> 500e05e2f8237f1cbc694e93414a4cc843c97dea
```
## Launching
Collect you samples in .fasta format in Samples dir and push in parent directory the command:
```
bash CovAnalyser.sh
<<<<<<< HEAD
=======

>>>>>>> 500e05e2f8237f1cbc694e93414a4cc843c97dea
```

## Output Files

| Extension | Description |
| --------- | ----------- |
<<<<<<< HEAD
| .csv | Tab-separated file of all features: Sample, Conclusion, Spike coverage, Spike mutations, Genetic Lines |
=======
| .csv | Tab-separated file of all features: Sample, Conclusion, Spike_coverage, Spike mutations, Genetic Lines |
>>>>>>> 500e05e2f8237f1cbc694e93414a4cc843c97dea
| .fa | File with sequence alignment of isolate and reference genome |

## Bugs

Submit problems or requests to the [Issue Tracker](https://github.com/artiman16/covanalyser/issues).

## Citation

# Licence