# CovAnalyser: SARS-CoV-2 genotyping program

## Introduction

CovAnalyser is a software tool to find mutations in SARS-CoV-2 Spike genes and determine genetic lineages of isolates.

## How does this work?

![Samples2](https://github.com/artiman16/covanalyser/assets/65973229/0d71ec8e-f0bc-489e-88c4-5f8f16237137)

## Installation (Linux)
```
git clone https://github.com/artiman16/covanalyser.git
cd covanalyser
conda env create -f environment.yml
```
## Launching
Collect you samples with .fasta format in Samples dir and push in parent directory the following commands:
```
conda activate covanalyser
bash CovAnalyser.sh
```

## Output Files

| Extension | Description |
| --------- | ----------- |
| .csv | Tab-separated file of all features: Sample, Conclusion, Spike_coverage, Spike mutations, Genetic Lines |
| .fa | File with sequence alignment of isolate and reference genome |

## Bugs

Submit problems or requests to the [Issue Tracker](https://github.com/artiman16/covanalyser/issues).

## Citation

Герасименко А. А., Водопьянов А. С., Писанов Р. В. Типирование штаммов SARS-COV-2 с помощью новой компьютерной программы «CovAnalyzer» // Сборник статей XXXVII международной научно-практической конференции «Российская наука в современном мире». М. – 2021. – С. 19-22.
 
# Licence
