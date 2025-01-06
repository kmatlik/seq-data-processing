# Code for processing sequencing data

This repository contains workflows for quality control, mapping, and visualisation of genomic datasets. Additional data processing steps for different types of data are listed below:
- ChIPseq: deduplication (picard) and peak calling (macs3)
- ATACseq: deduplication (picard) and peak calling (macs3)
- RNAseq: differential expression analysis based on counts over exons (DESeq2)
- nuclear RNAseq (FANSseq): differential expression analysis over exons as well as genes (DESeq2) 

The code is set up for R/RStudio installed on Windows, using WSL2 (Ubuntu) for packages that do not work on R or Windows (picard, macs3).

Some parts of the workflows are prepared with help from Rockefeller University Bioinformatics Resource Center training courses (https://rockefelleruniversity.github.io/).