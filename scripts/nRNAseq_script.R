# Load libraries
library(GenomicAlignments)
library(Rfastp)
library(Rsubread)
library(Rsamtools)
library(tidyverse)
library(magrittr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(DESeq2)
library(rtracklayer)

# Set up files & directories
read_1 <- c("R1_001.fastq.gz", "R1_002.fastq.gz", "R1_003.fastq.gz", "R1_004.fastq.gz")
outputdir <- "C:/Users/kvare/Documents/Seq/nRNAseq/Exp1"
samplelist <- c("sample1", "sample2", "sample3", "sample4")
samplegroups <- c("WT", "WT", "KO", "KO")
ref_genome <- "C:/Users/kvare/Documents/Seq/RefGenome/mm10_mainchrs"
id_to_symbol <- read_csv("C:/Users/kvare/Documents/Seq/mm10_id_to_symbol.csv")

# Quality control & trimming
for(i in 1:length(samplelist)){
  qc_report <- rfastp(
    read1 = read_1[i],
    outputFastq = paste0(outputdir, samplelist[i])
  )
}

# Mapping to reference genome
for(i in 1:length(samplelist)){
  
  # Align to reference genome
  mapping_report <- subjunc(index = ref_genome,
                            readfile1 = paste0(outputdir, samplelist[i], "_R1.fastq.gz"),
                            output_format = "BAM",
                            output_file = paste0(outputdir, samplelist[i], ".bam"), 
                            phredOffset = 33, 
                            nthreads = 8,
                            unique = TRUE,
                            useAnnotation = TRUE,
                            annot.inbuilt = "mm10")
  
  mapping_report %>% write.table(paste0(outputdir, samplelist[i], "_mapping_summary.csv"))
  
  # Sort and index BAM
  sortBam(file = paste0(outputdir, samplelist[i], ".bam"), 
          destination = paste0(outputdir, samplelist[i], "_sorted"))
  indexBam(paste0(outputdir, samplelist[i], "_sorted.bam"))
  
  # Remove unsorted bam files to save space
  file.remove(paste0(outputdir, samplelist[i], ".bam"))
}

# Counts over exons
mm10_exons <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
count_matrix <- matrix(data = 0, nrow = length(mm10_exons), ncol = length(samplelist))

for(i in 1:length(samplelist)){
  exon_counts <- summarizeOverlaps(features = mm10_exons,
                                   reads = paste0(outputdir, samplelist[i], "_sorted.bam"),
                                   ignore.strand = T)
  
  assay(exon_counts) %>% write.table(paste0(outputdir, samplelist[i], "_counts_exons.csv"))
  count_matrix[,i] <- assay(exon_counts)
}

rownames(count_matrix) <- names(mm10_exons)
colnames(count_matrix) <- samplelist
count_matrix %>% write.table(file = paste0(outputdir, "counts_exons.csv"))

# Counts over genebody
mm10_genes <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = "gene")
count_matrix_genebody <- matrix(data = 0, nrow = length(mm10_genes), ncol = length(samplelist))

for(i in 1:length(samplelist)){
  gene_counts <- summarizeOverlaps(features = mm10_genes,
                                   reads = paste0(outputdir, samplelist[i], "_sorted.bam"),
                                   ignore.strand = T)
  
  assay(gene_counts) %>% write.table(paste0(outputdir, samplelist[i], "_counts_genebody.csv"))
  count_matrix_genebody[,i] <- assay(gene_counts)
}

rownames(count_matrix_genebody) <- names(mm10_genes)
colnames(count_matrix_genebody) <- samplelist
count_matrix_genebody %>% write.table(file = paste0(outputdir, "counts_genebody.csv"))

# DESeq2 - exons
col_data <- data.frame(Group = samplegroups,
                       row.names = samplelist)

dds_exons <- DESeqDataSetFromMatrix(countData = count_matrix,
                                    colData = col_data,
                                    design = ~Group)
dds_exons <- DESeq(dds_exons)

res_exons <- results(dds_exons, contrast = c("Group", unique(samplegroups)[2], unique(samplegroups)[1]))
merge(id_to_symbol, res_exons %>% as.data.frame(), by.x = 1, by.y = 0) %>% write.table(file = paste0(outputdir, "DE_analysis_exons.csv"))
dds_exons %>% saveRDS(., file = paste0(outputdir, "dds_exons.RData"))

# DESeq2 - genebody
col_data <- data.frame(Group = samplegroups,
                       row.names = samplelist)

dds_genebody <- DESeqDataSetFromMatrix(countData = count_matrix_genebody,
                                       colData = col_data,
                                       design = ~Group)
dds_genebody <- DESeq(dds_genebody)

res_genebody <- results(dds_genebody, contrast = c("Group", unique(samplegroups)[2], unique(samplegroups)[1]))
merge(id_to_symbol, res_genebody %>% as.data.frame(), by.x = 1, by.y = 0) %>% write.table(file = paste0(outputdir, "DE_analysis_genebody.csv"))
dds_genebody %>% saveRDS(., file = paste0(outputdir, "dds_genebody.RData"))
