library(org.Hs.eg.db)
#library(Gviz)
#library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
# These are needed for EPplots.R
library(stringr)
library(dplyr)
library(GenomicRanges)
library(RColorBrewer)
library(beanplot)
#library(biomaRt)
library(parallel)

core_dir <- "src/core"
source(file.path(core_dir, "utils.R"))

path_file <- "src/core/paths.txt"
paths <- set_paths(path_file)
work_dir <- paths$work_dir
data_dir <- paths$data_dir
empP_dir <- paths$empP_dir
gwas_year <- paths$gwas_year
LDblock_path <- paths$LDblock_path
byGRBlist_name <- paths$byGRBlist_name

# Object that contains all GRBs, genes, enhancers and empP values for enh-transcript
byGRB_list <- readRDS(file.path(data_dir,"processed",byGRBlist_name))

loci_in_GRBs.df <- read.csv(file.path(data_dir,"processed","loci_in_GRBs.df.csv"), header = TRUE, stringsAsFactors = FALSE  )
# loci_in_GRBs.gr <- readRDS(file.path(data_dir,"processed","loci_in_GRBs.gr.rds"))
# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(data_dir,"inputs","hg19.canFam3.tsv"), sep = "\t", header = TRUE)  

TSSs.perTranscript.gr <- readRDS(file.path(data_dir, "inputs","TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(data_dir, "inputs","fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(data_dir, "inputs","promoter_cage_tpm_df_fixed_pooled.rds"))
# GRBS
grbs <- readRDS(file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))
# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(data_dir, "inputs","genes.hg19.ensembl.rds"))
gwas_associations <- read.csv(file.path(data_dir, "outputs","gwas_associations.csv"), header = TRUE, stringsAsFactors = FALSE)

# Create loci_in_GRB.gr from df
loci_in_GRBs.gr <- grangesFromDataframe(loci_in_GRBs.df, colName = LDblock, keep_metadata = TRUE)

# Filter loci_in_GRBs.gr to exclude those with NA in MAPPED_TRAIT
loci_in_GRBs.gr <- loci_in_GRBs.gr[!is.na(mcols(loci_in_GRBs.gr)$MAPPED_TRAIT)]

ep_pair_dir <- file.path("data/plots/expression_pairs")

# for testing one change value in bracket
# Parallel looping
traits <- unique(loci_in_GRBs.df$MAPPED_TRAIT)

# Define number of cores
num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

mclapply(traits, function(trait) {
    trait_dir <- str_replace_all(trait, " ", "_")
    trait_dir <- str_replace_all(trait_dir, "/", "_")
    if(!file.exists(file.path(ep_pair_dir,trait_dir))){
        dir.create(file.path(ep_pair_dir,trait_dir))
    }

    trait_grbs <- 
        filter(loci_in_GRBs.df, MAPPED_TRAIT == trait) %>% 
        pull(grb.coords) %>% 
        unique() %>%
        str_replace_all(pattern = ":|-", replacement = "_") 
    
    for(i_grb in trait_grbs){
        i_loci.gr <- loci_in_GRBs.gr[(mcols(loci_in_GRBs.gr)$MAPPED_TRAIT == trait)]
        plotEnhPromBeanPerTrait(byGRB_list, i=i_grb,
                                loci.gr = i_loci.gr,
                                neuronal=F, sign="permutation", 
                                dir=file.path(ep_pair_dir,trait_dir),
                                file_suffix= paste0(trait_dir,"_", i_grb, "_gene_expression_beanplot.png"),
                                enhancer_tpm.gr = fantom5.enhancers.tpm.gr,
                                promoter_tpm.df = promoter_cage_tpm_df,
                                verbose = TRUE)
    }
}, mc.cores = num_cores)
# # for (trait in unique(loci_in_GRBs.df$MAPPED_TRAIT)) {
  
#     # Generating folder in ep_pair_dir for trait
#     trait_dir <- str_replace_all(trait, " ", "_")
#     trait_dir <- str_replace_all(trait_dir, "/", "_")
#     if(!file.exists(file.path(ep_pair_dir,trait_dir))){
#         dir.create(file.path(ep_pair_dir,trait_dir))
#     }

#     # Grbs that are relevant for the trait
#     trait_grbs <- 
#         filter(loci_in_GRBs.df, MAPPED_TRAIT == trait) %>% 
#         pull(grb.coords) %>% 
#         unique() %>%
#         str_replace_all(pattern = ":|-", replacement = "_") 
    
#     for(i_grb in trait_grbs){
        
#         # loci_in_GRB for working trait
#         i_loci.gr <- loci_in_GRBs.gr[(mcols(loci_in_GRBs.gr)$MAPPED_TRAIT == trait)]
#         plotEnhPromBeanPerTrait(byGRB_list, i=i_grb,
#                                 # loci associated with trait
#                                 loci.gr = i_loci.gr,
#                                 neuronal=F, sign="permutation", 
#                                 dir=file.path(ep_pair_dir,trait_dir),
#                                 file_suffix= paste0(trait_dir,"_", i_grb, "_gene_expression_beanplot.png"),
#                                 enhancer_tpm.gr = fantom5.enhancers.tpm.gr,
#                                 promoter_tpm.df = promoter_cage_tpm_df,
#                                 verbose = FALSE)

#     }
# }

