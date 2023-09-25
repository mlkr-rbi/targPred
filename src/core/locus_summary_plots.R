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

source(file.path("utils.R"))

aux_dir <- "aux_files"


# Object that contains all GRBs, genes, enhancers and empP values for enh-transcript
byGRB_list <- readRDS(file.path("byGRB_list.rds"))

loci_in_GRBs.df <- readRDS(file.path("loci_in_GRBs.df.rds"))
loci_in_GRBs.gr <- readRDS(file.path("loci_in_GRBs.gr.rds"))
snps_in_GRBs.gr <- readRDS(file.path("snps_in_GRBs.gr.rds"))
# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(aux_dir,"hg19.canFam3.tsv"), sep = "\t", header = TRUE)  

TSSs.perTranscript.gr <- readRDS(file.path(aux_dir,"TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(aux_dir,"fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(aux_dir,"promoter_cage_tpm_df_fixed_pooled.rds"))
# GRBS
grbs <- readRDS(file.path("GRBs_with_targets.gr.rds"))
# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(aux_dir,"genes.hg19.ensembl.rds"))
gwas_associations <- readRDS(file.path(aux_dir,"gwas_associations.rds"))

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

col_vec <- c("darkgoldenrod1", "darkorange1", "red")
# plot the legend:
l <- c(expression("p">="0.01"), expression("p">="0.001"), expression("p"<="0.001"))
l_title <- "Empirical p-value"
plot_box_legend(l, col_vec, "legend_corrplot_empP.png", l_title, FALSE) # vertical legend
plot_box_legend(l, col_vec, "legend_corrplot_empP_horiz.png", l_title, TRUE) # horizontal legend


df <- byGRB_list[["chr2_29386917_30271314"]][[3]]
df["chr2:29435309-29435425","ENSG00000115295:ENST00000320081"]



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
  
  start_time <- Sys.time()
  
  for(i in 1:length(trait_grbs)){
  # For continuing grbs
    # loci_in_GRB for working trait
    i_loci.gr <- loci_in_GRBs.gr[(mcols(loci_in_GRBs.gr)$MAPPED_TRAIT == trait)]
    plots_prom_expression_beanplot_per_locus(byGRB_list, i=i,
                            neuronal=F, sign="permutation", 
                            dir=file.path(ep_pair_dir,trait_dir),
                            file_suffix= paste0(trait_dir, "_gene_expression_beanplot2.png"))
    
  }
  end_time <- Sys.time()
  print(start_time-end_time)
}

