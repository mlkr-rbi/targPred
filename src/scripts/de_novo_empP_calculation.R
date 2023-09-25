library(GenomicRanges)
library(plyranges)
library(ggpubr)
library(tidyverse)
library("effsize")


project_dir <- "/home/katarina/MyProjects/targPred"
inputs_dir <- "data/inputs"
processed_dir <- "data/processed"
outputs_dir <- "/data/outputs"

source(file.path(project_dir, "src/scripts/utils.R"))

# Max distance 6.4 Mb
max_dist <- 6400000

# Loading data
## These are pooled by cell/tissue type
fantom5.enhancers.tpm.gr <- readRDS(file.path(project_dir, processed_dir,"fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(project_dir, processed_dir,"promoter_cage_tpm_df_fixed_pooled.rds"))
transcripts <- readRDS(file.path(project_dir, processed_dir,  "TSSs.perTranscript.gr.rds"))

chr_vec <- paste0("chr", c(1:22)) %>%
  append("chrX")

EnhExpression <- function(fantom5.gr, enh.gr){
  return(apply(as.data.frame(mcols(subsetByOverlaps(fantom5.gr, enh.gr))), 2, sum))
}
EnhExpression2 <- function(enh_pos){
  enh.gr <- GRanges(enh_pos)
  return(apply(as.data.frame(mcols(subsetByOverlaps(fantom5.enhancers.tpm.gr, enh.gr))), 2, sum))
  
                                            
}

for(chr in chr_vec){
  print(paste("Processing chromosome", chr))
  
  # Subseting by chromosome
  enh_chr <- fantom5.enhancers.tpm.gr[seqnames(fantom5.enhancers.tpm.gr) == chr]
  transc_chr <- transcripts[seqnames(transcripts) == chr]
  
  # Preping some additional info
  transc_chr$gene <- paste0(transc_chr$ensembl_gene_id,":",transc_chr$ensembl_transcript_id)
  transc_chr$chr <- seqnames(transc_chr)
  transc_chr$transcript_position <- start(transc_chr)
  
  # Removing transcripts for which we do not have cage info
  transc_chr <- transc_chr[is.element(transc_chr$gene, rownames(promoter_cage_tpm_df))]
  
  # Converting to dfs for merge
  enh_chr.df <- as.data.frame(granges(enh_chr))
  
  # Large result
  empP_temp <- merge(enh_chr.df, mcols(transc_chr))
  
  # Filtering by max dist that will be tested
  empP_temp <- as.data.frame(empP_temp) %>%
    # end - start represnts middle of enhancers
    dplyr::filter(abs((end + start)/2 - transcript_position) <= max_dist) %>%
    mutate(enh_pos = paste0(seqnames, ":", start,"-",end,":",strand)) %>%
    dplyr::select(enh_pos, gene, transcript_position)
  
  empP_temp$empP <- NA
  # empP_temp <- grangesFromDataframe(empP_temp, colName = enh_pos)
  
  # Vectorizing
  enh_expression.ls <- vector(mode = "list", length = length(unique(empP_temp$enh_pos)))
  names(enh_expression.ls) <- unique(empP_temp$enh_pos)[1:5]
  
  enh_expression.ls <- lapply(unique(empP_temp$enh_pos)[1:5], EnhExpression2)
  names(enh_expression.ls) <- unique(empP_temp$enh_pos)[1:5]
  
  ApplyExpressionPair <- function(index1, index2, enh.ls, promoter.df){
    print(index1)
    print( enh.ls[[index1]])
    
    expression_pair <- data.frame(enh_expression = enh.ls[[index1]],
                                  prom = promoter.df[index2])
    return(empP_diff_in_median_sampling(expression_pair, iter = 100))
  }
  
  temp <- as.data.frame(mcols(empP_temp))[1:5,] %>%
    mutate(empP = ApplyExpressionPair(enh_pos, gene, enh_expression.ls, promoter_cage_tpm_df))
  
  apply(as.data.frame(mcols(empP_temp)), 1, enh_expression.ls)
  
  first_na <- which(is.na(enh_temp$empP), arr.ind = TRUE)[1]
  
  for( i in first_na:length(empP_temp)){
    
      enh_expression <- EnhExpression(fantom5.gr = fantom5.enhancers.tpm.gr,
                                      enh.gr = empP_temp[i])
      
      empP_temp[i]$empP <- 
        empP_diff_in_median_sampling(data.frame(
          enh_expression = enh_expression,
          prom = promoter_cage_tpm_df[rownames(promoter_cage_tpm_df) == empP_temp[i]$gene]))
  }
  if(i %% 1000 == 0){
    writeGrangesToText(empP_temp, is_empP_file = TRUE, dir = file.path(project_dir, processed_dir, paste0(chr, "empP_files_denovo/","empP_CAGEen_denovo_hg19Ens_6.4Mb.rds")))    
  }
  
    
}
