library(Rfastp)
library(Rsubread)
library(Rsamtools)
library(tidyverse)
library(magrittr)
library(rtracklayer)

#Next:
#ATAC and ChIP counts over peaks
#ATAC consensus peaks
#ATAC DE peaks

#Test arguments
args <- c("C:/Users/kvare/Documents/Seq/sample_info.xlsx", "C:/Users/kvare/Documents/Seq/ATACseq/Exp1/", "ATAC") #path to sampleinfo, outputdir

#args <- commandArgs(TRUE)

#args[1] = path to sample sheet (Excel file)
#args[2] = output directory
#args[3] = one of the following: "ATAC", "ChIP", "nRNA", "RNA"

if (length(args)<3) {
  
  #If some arguments are not supplied
  print("Arguments missing")
  
}else{

# Import sample information  
  data_table <- readxl::read_xlsx(path = args[1]) 
  read_1 <- data_table$Read1
  read_2 <- data_table$Read2
  samplelist <- data_table$SampleName
  inputlist <- data_table$InputName
  samplegroups <- data_table$SampleGroup

  outputdir <- args[2]
  outputdir_wsl <- paste0("/mnt/c", strsplit(args[2], split = ":")[[1]][2]) #Changing to wsl compatible path. This is required to run picard and macs3 through wsl. Note that the path starts with /mnt/c/ instead of C:/ 
  ref_genome <- "C:/Users/kvare/Documents/Seq/RefGenome/mm10_mainchrs"
  id_to_symbol <- read_csv("C:/Users/kvare/Documents/Seq/mm10_id_to_symbol.csv")

# Quality control & trimming
  for(i in 1:length(samplelist)){
    qc_report <- rfastp(
      read1 = read_1[i],
      read2 = ifelse(args[3]=="ATAC", read_2[i], ""),
      outputFastq = paste0(outputdir, samplelist[i])
    )
  }

# Mapping to reference genome
   for(i in 1:length(samplelist)){
    
    # Align to reference genome
    if(args[3]=="ATAC"){
      
      mapping_report <- align(index = ref_genome,
                              readfile1 = paste0(outputdir, samplelist[i], "_R1.fastq.gz"),
                              readfile2 = paste0(outputdir, samplelist[i], "_R2.fastq.gz"),
                              output_format = "BAM",
                              output_file = paste0(outputdir, samplelist[i], ".bam"), 
                              type = "dna",
                              phredOffset = 33, 
                              nthreads = 4,
                              unique = TRUE,
                              maxFragLength = 2000)
    
      }else if(args[3]=="ChIP"){
      
        mapping_report <- align(index = ref_genome,
                              readfile1 = paste0(outputdir, samplelist[i], "_R1.fastq.gz"),
                              output_format = "BAM",
                              output_file = paste0(outputdir, samplelist[i], ".bam"), 
                              type = "dna",
                              phredOffset = 33, 
                              nthreads = 4,
                              unique = TRUE)
    
        }else if(args[3]=="RNA" | args[3]=="nRNA"){ #RNA or nRNA
          
          mapping_report <- subjunc(index = ref_genome,
                                    readfile1 = paste0(outputdir, samplelist[i], "_R1.fastq.gz"),
                                    output_format = "BAM",
                                    output_file = paste0(outputdir, samplelist[i], ".bam"), 
                                    phredOffset = 33, 
                                    nthreads = 4,
                                    unique = TRUE,
                                    useAnnotation = TRUE,
                                    annot.inbuilt = "mm10")
      
    }
    
    #Write a report
    mapping_report %>% write.table(paste0(outputdir, samplelist[i], "_mapping_summary.csv"))
  
    # Sort BAM
    sortBam(file = paste0(outputdir, samplelist[i], ".bam"), 
            destination = paste0(outputdir, samplelist[i], "_sorted"))
  
    # Remove unsorted bam files to save space
    file.remove(paste0(outputdir, samplelist[i], ".bam"))
  }

  
  
# Run picard
  if(args[3]=="ATAC" | args[3]=="ChIP"){
    
    # Run picard AddOrReplaceReadGroups to add read groups
    for(i in 1:length(samplelist)){
      picard_command <- paste0("wsl java -jar /home/kart/picard/build/libs/picard.jar AddOrReplaceReadGroups -I ",
                               outputdir_wsl, samplelist[i], "_sorted.bam",
                               " -O ",
                               outputdir_wsl, samplelist[i], "_withRG.bam",
                               " -ID 1 -LB lib1 -PL illumina -PU unit1 -SM 1")
      system(picard_command)
    }
    
    # Remove extra bam files to save space
    for(i in 1:length(samplelist)){
      file.remove(paste0(outputdir, samplelist[i], "_sorted.bam")) #Delete sorted bam files without read groups
    }
    
    # Run picard MarkDuplicates. Generate bam files that are sorted and deduplicated.
    for(i in 1:length(samplelist)){
      picard_command <- paste0("wsl java -jar /home/kart/picard/build/libs/picard.jar MarkDuplicates -I ",
                               outputdir_wsl, samplelist[i], "_withRG.bam",
                               " -O ",
                               outputdir_wsl, samplelist[i], "_sorted.bam",
                               " -M ",
                               outputdir_wsl, samplelist[i], "_metrix.csv, --REMOVE_DUPLICATES true")
      system(macs_command)
    }
    
    # Remove extra bam files to save space
    for(i in 1:length(samplelist)){
      file.remove(paste0(outputdir, samplelist[i], "_withRG.bam"))
    }
  }
  
# Index bams
  for(i in 1:length(samplelist)){
    indexBam(paste0(outputdir, samplelist[i], "_sorted.bam"))
  }

# Create bigwigs
  for(i in 1:length(samplelist)){
    # Mapping summary over chromosomes
    mapped_stats <- idxstatsBam(paste0(outputdir, samplelist[i], "_sorted.bam"))
    mapped_stats %>% write_csv(paste0(outputdir, samplelist[i], "_mapped_stats.csv"))
    
    # Create bigwig files
    read_coverage <- coverage(paste0(outputdir, samplelist[i], "_sorted.bam"),
                              weight = (10^6)/sum(mapped_stats[,"mapped"]))
    export.bw(read_coverage, paste0(outputdir, samplelist[i], ".bw"))
  }

# Peak calling with macs3
  if(args[3]=="ATAC"){
    for(i in 1:length(samplelist)){
      macs_command <- paste0(c("wsl /home/kart/miniconda3/bin/macs3 callpeak -t "), 
                             outputdir_wsl, samplelist[i], "_sorted.bam",
                             " -n ",
                             samplelist[i],
                             " -f BAMPE --outdir ",
                             outputdir_wsl,
                             " -g mm")
      system(macs_command)
    }
    
  }else if(args[3]=="ChIP"){
    for(i in 1:length(samplelist)){
      macs_command <- paste0(c("wsl /home/kart/miniconda3/bin/macs3 callpeak -t "), 
                             outputdir_wsl, samplelist[i], "_sorted.bam",
                             " -c ",
                             outputdir_wsl, inputlist[i], "_sorted.bam",
                             " -n ",
                             samplelist[i],
                             " -f BAM --outdir ",
                             outputdir_wsl,
                             " -g mm")
      system(macs_command)
    }
  }
  
# Counts over features
  if(args[3]=="nRNA"){
    
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
  
  }else if(args[3]=="RNA"){
    
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
  
  }else if(args[3]=="ATAC"){
    
    # import consensus peaks and count over consensus peaks
  }
  
# Differential expression analysis
  if(args[3]=="nRNA"){
    
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
  
  }else if(args[3]=="RNA"){
      
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
    
  }else if(args[3]=="ATAC"){
      
    #Add ATAC DE analysis
  }
  
  

}
