library(org.Hs.eg.db)
library(Gviz)
#library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
# These are needed for EPplots.R
library(stringr)
library(dplyr)
library(GenomicRanges)
library(RColorBrewer)
library(beanplot)
library(biomaRt)

path_file <- "paths.txt"
path_vector <- readr::read_lines("paths.txt")
directory_type <- ""
core_dir <- ""
work_dir <- ""
data_dir <- ""
gwas_year <- ""
byGRBlist_name <- ""

for (idx in 1:length(path_vector)) {
    line_vector <- strsplit(path_vector[idx],"")
    if (line_vector[[1]][1] == "$") {
        directory_type <- paste(line_vector[[1]][-1], collapse = "")
        next
    } else if (line_vector[[1]][1] == "#") {
        # we can comment a line with symbol #
        next
    } else {
        if (directory_type == "core") core_dir <- path_vector[idx]
        if (directory_type == "work") work_dir <- path_vector[idx]
        if (directory_type == "data") data_dir <- path_vector[idx]
        if (directory_type == "year") gwas_year <- path_vector[idx]
        if (directory_type == "byGRB_name") byGRBlist_name <- path_vector[idx]
    }
}

if (core_dir == "") {
    print("ERROR: core_dir not defined in paths.txt!")
    quit(status=1)
} else if (work_dir == "") {
    print("ERROR: work_dir not defined in paths.txt!")
    quit(status=1)
} else if (data_dir == "") {
    print("ERROR: data_dir not defined in paths.txt!")
    quit(status=1)
} else if (gwas_year == "") {
    print("ERROR: gwas_year not defined in paths.txt!")
    quit(status=1)
}
if (byGRBlist_name == "") {
    print("byGRB_list name is default byGRB_list.rds")
    byGRBlist_name <- "byGRB_list.rds"
}

source(file.path(core_dir, "utils.R"))

# Object that contains all GRBs, genes, enhancers and empP values for enh-transcript
byGRB_list <- readRDS(file.path(data_dir,"outputs",byGRBlist_name))

loci_in_GRBs.df <- readRDS(file.path(data_dir,"processed","loci_in_GRBs.df.rds"))
loci_in_GRBs.gr <- readRDS(file.path(data_dir,"processed","loci_in_GRBs.gr.rds"))
# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(data_dir,"inputs","hg19.canFam3.tsv"), sep = "\t", header = TRUE)  

TSSs.perTranscript.gr <- readRDS(file.path(data_dir, "processed","TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(data_dir, "processed","fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(data_dir, "processed","promoter_cage_tpm_df_fixed_pooled.rds"))
# GRBS
grbs <- readRDS(file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))
# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(data_dir, "processed","genes.hg19.ensembl.rds"))
gwas_associations <- readRDS(file.path(data_dir, "inputs", paste0("gwas_catalog_",gwas_year),"gwas_associations.rds"))

ep_pair_dir <- file.path("plots/expression_pairs")

# for testing one change value in bracket
for (trait in unique(loci_in_GRBs.df$MAPPED_TRAIT)) {
  
    # Generating folder in ep_pair_dir for trait
    trait_dir <- str_replace_all(trait, " ", "_")
    trait_dir <- str_replace_all(trait_dir, "/", "_")
    if(!file.exists(file.path(ep_pair_dir,trait_dir))){
        dir.create(file.path(ep_pair_dir,trait_dir))
    }

    # Grbs that are relevant for the trait
    trait_grbs <- 
        filter(loci_in_GRBs.df, MAPPED_TRAIT == trait) %>% 
        pull(grb.coords) %>% 
        unique() %>%
        str_replace_all(pattern = ":|-", replacement = "_") 
    
    for(i_grb in trait_grbs){
        
        # loci_in_GRB for working trait
        i_loci.gr <- loci_in_GRBs.gr[(mcols(loci_in_GRBs.gr)$MAPPED_TRAIT == trait)]
        plotEnhPromBeanPerTrait(byGRB_list, i=i_grb,
                                # loci associated with trait
                                loci.gr = i_loci.gr,
                                neuronal=F, sign="permutation", 
                                dir=file.path(ep_pair_dir,trait_dir),
                                file_suffix= paste0(trait_dir,"_", i_grb, "_gene_expression_beanplot.png"),
                                enhancer_tpm.gr = fantom5.enhancers.tpm.gr,
                                promoter_tpm.df = promoter_cage_tpm_df)

    }
}

