# Code for processing sequencing data

## Intro
This repository contains workflows for quality control (Rfastp), mapping (Rsubread), and visualisation of genomic datasets. Additional data processing steps for different types of data are listed below:
- ATACseq (paired-end reads): deduplication (picard) and peak calling (macs3)
- ChIPseq (single-end): deduplication (picard) and peak calling (macs3)
- RNAseq: differential expression analysis based on counts over exons (DESeq2)
- nuclear RNAseq (FANSseq): differential expression analysis over exons as well as genes (DESeq2) 

The code is set up for R/RStudio installed on Windows, using WSL2 (Ubuntu) for packages that do not work on R or Windows (picard, macs3).

Part of the workflows is prepared with help from Rockefeller University Bioinformatics Resource Center training courses (https://rockefelleruniversity.github.io/).

## Sample information
Sample information should be provided as an Excel file with the following column names:
SampleName = list of samples, including ChIPseq input samples  
SampleGroup = experimental group, required for DESeq2 (ATACseq, RNAseq, and nRNAseq)  
Read1 = path to Read 1  
Read2 = path to Read 2, required for ATACseq  
InputName = sample to be used as input, required for ChIPseq peak calling with MACS3  

Example:
 | SampleName | SampleGroup | Read1           | Read2           | InputName |
 | ---------- | ----------- | --------------- | --------------- | --------- |
 | sample1    | group1      | R1_001.fastq.gz | R2_001.fastq.gz | sample4   |
 | sample2    | group1      | R1_002.fastq.gz | R2_002.fastq.gz | sample4   |
 | sample3    | group2      | R1_003.fastq.gz | R2_003.fastq.gz | sample4   |
 | sample4    | group2      | R1_004.fastq.gz | R2_004.fastq.gz | sample4   |

## Running the script
Code:
```
Rscript --vanilla myscript.R arg1 arg2 arg3
```
- arg1 = absolute path to sample sheet (xlsx file) 
- arg2 = absolute path to output directory 
- arg3 = type of sequencing data; one of the following: "ATAC", "ChIP", "RNA", "nRNA" 

Example: 
```
Rscript --vanilla myscript.R "C:/Users/kvare/Documents/Sequencing/sample_info.xlsx" "C:/Users/kvare/Documents/Sequencing/Exp1/" "ATAC"
```
