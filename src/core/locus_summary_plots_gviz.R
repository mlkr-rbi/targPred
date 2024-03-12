library(org.Hs.eg.db)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
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
snps_in_GRBs.gr <- readRDS(file.path(data_dir,"processed","snps_in_GRBs.gr.rds"))
# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(data_dir,"inputs","hg19.canFam3.tsv"), sep = "\t", header = TRUE)  
# Building target gene genomic ranges
# enhancers and promoters are pooled data
TSSs.perTranscript.gr <- readRDS(file.path(data_dir, "processed","TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(data_dir, "processed","fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(data_dir, "processed","promoter_cage_tpm_df_fixed_pooled.rds"))
# GRBS
grbs <- readRDS(file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))
# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(data_dir, "processed","genes.hg19.ensembl.rds"))
gwas_associations <- readRDS(file.path(data_dir, "inputs", paste0("gwas_catalog_",gwas_year),"gwas_associations.rds"))

col_vec <- c("darkgoldenrod1", "darkorange1", "red")

gviz_dir <- file.path( "plots/gviz_plots")

# Adding a column to grbs that correspond to targets
# actual plotting, saved in gviz_dir
loci_in_GRBs.gr$name <- loci_in_GRBs.gr$SNPS # need to have a $name column for gviz plotting

# preparing mart to extract gene objects from
mart <- useMart(dataset="hsapiens_gene_ensembl", biomart="ENSEMBL_MART_ENSEMBL", host="https://feb2014.archive.ensembl.org")


for(trait in unique(loci_in_GRBs.df$MAPPED_TRAIT)){
    trait_dir <- file.path(gviz_dir,str_replace_all(trait," |/","_"))
    # Creates a dir if it does not exist
    if(!file.exists(trait_dir)){dir.create(trait_dir)}
 
    # to be done on each locus
    message(trait)
    # Including only loci in trait
    loci_in_GRB_trait <- loci_in_GRBs.gr[loci_in_GRBs.gr$MAPPED_TRAIT == trait,]
    # Remoivng duplicated SNPS per trait/locus
    loci_in_GRB_trait <- loci_in_GRB_trait[!duplicated(loci_in_GRB_trait$SNPS),]
    
    grbs_to_loop <-
        which( names(byGRB_list) %in% str_replace_all(loci_in_GRB_trait$grb.coords, pattern = ":|-", replacement = "_"))
    
    #for (i in (1:length(byGRB_list))) {
    for (i in grbs_to_loop) {
        print(paste("Plotting Gviz for", names(byGRB_list)[i], sep=" "))
        tryCatch(locusGviz(i, 
                           dir=paste0(trait_dir,"/",str_replace_all(trait," |/","_"),"_"),
                           labels=T, 
                           byGRB_list=byGRB_list, 
                           # Selecting only loci with the trait
                           loci.gr= loci_in_GRB_trait, 
                           # Selecting only SNPs with the trait
                           snps.gr=snps_in_GRBs.gr[snps_in_GRBs.gr$SNPS %in% loci_in_GRB_trait$SNPS,],
                           grbs=grbs, 
                           mart=mart,
                           trait = trait),
                 error=function(e){cat("ERROR :",conditionMessage(e), "\n", trait, i)})
    }
}

