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
for (trait in unique(loci_in_GRBs.df$MAPPED_TRAIT)[1]) {
  
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
  #for(i_grb in trait_grbs[match(i_grb, trait_grbs):length(trait_grbs)]){
      
    # loci_in_GRB for working trait
    i_loci.gr <- loci_in_GRBs.gr[(mcols(loci_in_GRBs.gr)$MAPPED_TRAIT == trait)]
    plots_prom_expression_beanplot_per_locus(byGRB_list, i=i,
                            neuronal=F, sign="permutation", 
                            dir=file.path(ep_pair_dir,trait_dir),
                            file_suffix= paste0(trait_dir,"_", "_gene_expression_beanplot2.png"))
#    break
    
  }
  end_time <- Sys.time()
  print(start_time-end_time)
#  break
}

#comment out gviz as it is done in another script
if (FALSE){
gviz_dir <- file.path( "plots/gviz_plots")


# Adding a column to grbs that correspond to targets
grbs$targets <- grbs$geneName
# actual plotting, saved in gviz_dir
loci_in_GRBs.gr$name <- loci_in_GRBs.gr$SNPS # need to have a $name column for gviz plotting
# pretty sure this is not important. transfered this into GWAS_SNPs_inGRB
# genome(loci_in_GRBs.gr) <- genome(BSgenome.Hsapiens.UCSC.hg19)
# seqlevels(loci_in_GRBs.gr) <- seqlevels(BSgenome.Hsapiens.UCSC.hg19)
# seqlengths(loci_in_GRBs.gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg19)
# gen <- genome(loci_in_GRBs.gr)
# preparing mart to extract gene objects from
mart <- useMart(dataset="hsapiens_gene_ensembl", biomart="ENSEMBL_MART_ENSEMBL", host="https://feb2014.archive.ensembl.org")


#gging loci to be treated separately as these contain some exceptions:
# TBR1 locus [10] -> merge 2 GRBs...
#sink("Gviz_Rconsoleoutput.txt")
for(trait in unique(loci_in_GRBs.df$MAPPED_TRAIT)){
  trait_dir <- file.path(gviz_dir,str_replace_all(trait," |/","_"))
  # Creates a dir if it does not exist
  if(!file.exists(trait_dir)){dir.create(trait_dir)}
  #sink(paste(gviz_dir,"/",str_replace_all(trait," |/","_"),"/","Gviz_Rconsoleoutput.txt", sep = ""))
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
                snps.gr=snps_in_GRBs.gr[snps_in_GRBs.gr$SNPS %in%
                                          loci_in_GRB_trait$SNPS,],
                grbs=grbs, 
                mart=mart,
                trait = trait),
              error=function(e){cat("ERROR :",conditionMessage(e), "\n",
                                    trait, i)})
  }
  #sink()
}


Purple <- "#A74F7C" 
Grey   <- "#7D949B"
print(snps_in_GRBs.gr[snps_in_GRBs.gr$SNPS %in% c("rs1891058", "rs147340715")])
current_SNPs <- 
#  filter(snps_in_GRBs.gr, is.element(SNPS, c("rs1891058", "rs147340715")))
  snps_in_GRBs.gr[snps_in_GRBs.gr$SNPS %in% c("rs1891058", "rs147340715")]
current_GRB <- subsetByOverlaps(grbs, current_SNPs)
strack <- AnnotationTrack(current_SNPs, name="GWAS SNPs", 
                          stacking="full",
                          group=current_SNPs$SNPS, groupAnnotation="group",
                          fill=Purple, col=Purple)
grbtrack <- AnnotationTrack(current_GRB, name="GRB human:dog", 
                            id=current_GRB$targets, fontcolor.feature="black", 
                            fill=Grey, col=Grey)
displayPars(grbtrack) <- list(featureAnnotation="id")
plotTracks(list(grbtrack, strack),
           sizes = c(1,1),
           background.title="azure3", fontcolor.title="black")
}
