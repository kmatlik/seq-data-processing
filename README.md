# Code for processing sequencing data

## Intro
This repository contains workflows for quality control (Rfastp), mapping (Rsubread), and visualisation of genomic datasets. Additional data processing steps for different types of data are listed below:
- ChIPseq: deduplication (picard) and peak calling (macs3)
- ATACseq: deduplication (picard) and peak calling (macs3)
- RNAseq: differential expression analysis based on counts over exons (DESeq2)
- nuclear RNAseq (FANSseq): differential expression analysis over exons as well as genes (DESeq2) 

The code is set up for R/RStudio installed on Windows, using WSL2 (Ubuntu) for packages that do not work on R or Windows (picard, macs3).

Part of the workflows is prepared with help from Rockefeller University Bioinformatics Resource Center training courses (https://rockefelleruniversity.github.io/).

## Running the script
Code:
```
Rscript --vanilla myscript.R arg1 arg2 arg3
```
arg1 = absolute path to sample sheet (xlsx file) 
arg2 = absolute path to output directory 
arg3 = type of sequencing data; one of the following: "ATAC", "ChIP", "RNA", "nRNA" 

Example: 
```
Rscript --vanilla myscript.R "C:/Users/kvare/Documents/Sequencing/sample_info.xlsx" "C:/Users/kvare/Documents/Sequencing/Exp1/" "ATAC"
```
