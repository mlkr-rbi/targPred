#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager",repos = "http://cran.us.r-project.org")
#if (!require("tidyverse", quietly = TRUE))
#    install.packages("tidyverse",repos = "http://cran.us.r-project.org")
#if (!require("acepack", quietly = TRUE))
#    install.packages("acepack",repos = "http://cran.us.r-project.org")

#BiocManager::install("rtracklayer")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("Gviz")
#BiocManager::install("biovizBase")
#BiocManager::install("plyranges")
#BiocManager::install("BSgenome")
#BiocManager::install("S4Vectors")
#BiocManager::install("IRanges")
#BiocManager::install("GenomicRanges")


library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(rols)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
# library(purrr)
# library(curl)
# library(RColorBrewer)
# library(org.Hs.eg.db)
# library(Gviz)
# library(BSgenome.Hsapiens.UCSC.hg19)

# library(httr)
# library(jsonlite)
# library(xml2)

core_dir <- "/home/katarina/MyProjects/targPred/src/core"
data_dir <- "/home/katarina/MyProjects/targPred/data"

source(file.path(core_dir, "utils.R"))



# Loading data
aux_dir <- "../aux_files"

# Output from gwas_associations.py
 gwas <- read.csv(file.path(data_dir, "outputs", "gwas_associations_LDblock.csv"))
   
 
# Main GRB dataset
cne_density_grb <- 
  import.bed(file.path(data_dir, "inputs","hg19.canFam3.98.50.processed.threshold.grbs.bed")) %>%
  # Reducing since some GRBs overlap- Look at CNE_density_overlapping_GRBs.txt for examples
  GenomicRanges::reduce()

# Modifing because granges have 1 base and bed files are 0 based
start(cne_density_grb) = start(cne_density_grb) -1

# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(data_dir, "inputs","hg19.canFam3.tsv"), sep = "\t", header = TRUE)


# Helper functions
grangesFromDataframe <- function(data, colName, keep_metadata=TRUE){
  # Description: make granges from dataframe
  # which has a column of this format:
  # chr4:11394774-11431389
  data_r <-
    data %>%
    mutate(start=as.character(str_split(!!rlang::enquo(colName), pattern = ":|-", simplify = TRUE)[,2]),
           end=as.character(str_split(!!rlang::enquo(colName), pattern = ":|-", simplify = TRUE)[,3]),
           CHR= as.character(str_split(!!rlang::enquo(colName), pattern = ":|-", simplify = TRUE)[,1]))
  
  
  # Granges object is needed to detect SNP loci that overlap 2 GRBs
  data_r <-
    makeGRangesFromDataFrame(data_r,
                            seqnames.field ="chr",
                            start.field = "start",
                            end.field = "end",
                            keep.extra.columns = keep_metadata)
  return(data_r)
}

# Preprocessing gwas data

non_gene_names <- c("",NA,"intergenic","Intergenic", "NR")
gwas <- 
  gwas %>%
  select(CHR= "CHR_ID", SNPS, LDblock, STUDY.ACCESSION, REPORTED.GENE.S., 
         MAPPED_GENE, P.VALUE, MAPPED_TRAIT, MAPPED_TRAIT_URI) %>%
  filter(!is.na(CHR))
gwas <- gwas %>%
  # Creating a temporary column that contains all genes comma sep
  mutate(temp = case_when(is.element(REPORTED.GENE.S., non_gene_names) ~ MAPPED_GENE,
                                is.element(MAPPED_GENE, non_gene_names) ~ REPORTED.GENE.S.,
                                TRUE ~ str_c(REPORTED.GENE.S., MAPPED_GENE, sep = ", "))) %>%
  # Removing repeated genes (eg. SKIV2L, SKIV2L, SKIV2L, SKIV2L, SKIV2L, SKIV2L)
  mutate(GWAS_GENES = map_chr(str_split(temp, ", "), ~str_c(unique(.x),collapse = ", "))) %>%
  dplyr::select(-c(temp, REPORTED.GENE.S., MAPPED_GENE))

## gwas_gr stores SNP LD blocks in GenomicRanges object
gwas_gr <- grangesFromDataframe(gwas, LDblock)


# Building target gene GRanges
target_genes_gr <- grangesFromDataframe(target_genes,
                                        colName =ucscCoordinates.hg19)

# Adding target genes as metadata to grb
cne_density_grb <- 
  plyranges::join_overlap_left(cne_density_grb, target_genes_gr) %>%
  # Using dplyr to merge multiple targets in GRBs
  as.data.frame() %>%
  dplyr::group_by(start, end) %>%
  # Concatenating all these cols from metadata
  mutate(ENSEMBLGeneID = paste0(ENSEMBLGeneID, collapse = ", "),
         geneName = paste0(geneName, collapse = ", "),
         originalSupport = paste0(originalSupport, collapse = ", "),
         normalisedSupport = paste0(normalisedSupport, collapse = ", "),
         ucscCoordinates.hg19 = paste0(ucscCoordinates.hg19, collapse = ", ")) %>%
  unique()%>%
  ungroup() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Overlaping gwas loci with GRBs
loci_in_GRBs.gr <- subsetByOverlaps(gwas_gr, cne_density_grb)
ovlp <- findOverlaps(loci_in_GRBs.gr, cne_density_grb)

# Adding GRB coordinates to each locus as character
loci_in_GRBs.gr$grb.coords[queryHits(ovlp)] <-
  as.character(cne_density_grb[subjectHits(ovlp)])

# Adding GRB target genes 
loci_in_GRBs.gr$GRB_GENES[queryHits(ovlp)] <- 
  cne_density_grb[subjectHits(ovlp)]$geneName

# Saving loci in GRB as bed file for visualizations in Genome browser
export.bed(loci_in_GRBs.gr, con = file.path(data_dir, "outputs", "loci_in_GRBs.bed"))


loci_in_GRBs.df <- mcols(loci_in_GRBs.gr) %>%
  as.data.frame()

# Separating Traits
## Trait ontology contains all GWAS trait EFO IDs extracted from MAPPED_TRAIT_URI
trait_ontology <- dplyr::select(loci_in_GRBs.df, MAPPED_TRAIT_URI) %>%
  unique()
rownames(trait_ontology) <- NULL
trait_ontology <- 
  trait_ontology %>%
  # Separating multiple traits into rows
  tidyr::separate_rows(MAPPED_TRAIT_URI, sep = ",") %>%
  mutate(MAPPED_TRAIT_URI = stringr::str_trim(MAPPED_TRAIT_URI, side = "both")) %>%
  unique() %>%
  # Splitting the URI eg. http://www.ebi.ac.uk/efo/EFO_0007978 to select only the ID EFO_0007978
  mutate(EFO_ID = str_split(MAPPED_TRAIT_URI, pattern = "/", simplify = TRUE)[,5]) %>%
  filter(!is.na(MAPPED_TRAIT_URI)) %>%
  filter(MAPPED_TRAIT_URI != "")


# Downloaded efo terms
efo_terms <- readRDS(file.path(data_dir, "inputs","EFO_terms.rds"))


## Relabeling Traits
# Using trait_ontology dataframe to do this once for each term
trait_ontology$label <- NA

for(id in 1:length(trait_ontology$EFO_ID)){
  
  term_id <- str_replace(trait_ontology$EFO_ID[id], "_",":")
  
  # efo_terms is an ontology object which needs the termLabel method to access label
  tryCatch(trait_ontology$label[id] <- rols::termLabel(efo_terms[term_id]),
           error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}



loci_in_GRBs.df <-
  loci_in_GRBs.df %>%
  # Multiple URIs that correspond to one trait are sep with ","
  separate_rows(MAPPED_TRAIT_URI, sep = ",") %>%
  mutate(MAPPED_TRAIT_URI = stringr::str_trim(MAPPED_TRAIT_URI, side = "both")) %>%
  # Splitting URI slashes since the unique ID is stored in the end
  mutate(MAPPED_TRAIT_ID = str_split(MAPPED_TRAIT_URI, pattern = "/", simplify =TRUE)[,5]) %>%
  left_join(dplyr::select(trait_ontology, MAPPED_TRAIT_URI, label)) %>%
  filter(MAPPED_TRAIT_URI != "")



# Removing the obsolete_ addition from some trait labels
loci_in_GRBs.df <- 
  loci_in_GRBs.df %>%
  mutate(MAPPED_TRAIT = ifelse(str_detect(MAPPED_TRAIT,"obsolete_"),
                               str_remove(MAPPED_TRAIT, "obsolete_"), 
                               MAPPED_TRAIT)) %>%
  # Removing NA MAPPED_TRAIT
  filter(!is.na(MAPPED_TRAIT))
  



# SNPs in GRBs granges

# Coercing to data frame because 
# it makes a GPos object by default which makes bugs afterwards
# snps_in_GRBs.gr <- 
#   snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, 
#            unique(loci_in_GRBs.df$SNPS),
#            ifnotfound = "drop") %>%
#   as.data.frame() %>%
#   dplyr::select(seqnames, pos, SNPS = "RefSNP_id")
# # Generating Granges 
# snps_in_GRBs.gr <- 
#   makeGRangesFromDataFrame(snps_in_GRBs.gr,
#                            seqnames.field = "seqnames",
#                            start.field = "pos",
#                            end.field ="pos",
#                            keep.extra.columns = TRUE)


# Concetinating  to the final object
## Changing type of seqnames from 1,2,3,... to chr1, chr2, chr3,...
# seqlevelsStyle(snps_in_GRBs.gr) <- "UCSC"
# 
# # ?Maybe combine with loci_in_GRBs.df for extra metadata
# saveRDS(snps_in_GRBs.gr, file = file.path("snps_in_GRBs.gr.rds"))

####### byGRB_list ##########

# Loading data

# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(data_dir, "inputs","hg19.canFam3.tsv"), sep = "\t", header = TRUE)  

# Building target gene genomic ranges
target_genes_gr <- grangesFromDataframe(target_genes, ucscCoordinates.hg19)

# Sent by Anja by mail (from fantom5 data)
TSSs.perTranscript.gr <- readRDS(file.path(data_dir, "processed","TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(data_dir, "processed","fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(data_dir, "processed","promoter_cage_tpm_df_fixed_pooled.rds"))

grbs <- cne_density_grb

# Adding target genes to GRBs 
grbs$targets <- ""
for (i in 1:length(grbs)) {
  t <- sort(subsetByOverlaps(target_genes_gr, grbs[i]))
  if(length(t)>0) {
    grbs$targets[i] <- paste(t$geneName, collapse=",")
  }
}
saveRDS(grbs, file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))

# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(data_dir, "processed","genes.hg19.ensembl.rds"))
gwas_associations <- readRDS(file.path(data_dir, "inputs/gwas_catalog_2021","gwas_associations.rds"))


# makes sure to work on the GRB level, even when one GRB covers multiple loci
grb_subset <- subsetByOverlaps(grbs, loci_in_GRBs.gr)

# add grb targets to loci in grbs
grb_subset.df <- data.frame(chr=seqnames(grb_subset),
                            start=start(grb_subset),
                            end=end(grb_subset),
                            targets=grb_subset$targets,
                            stringsAsFactors=F)

grb_subset.df <-
  grb_subset.df %>%
  # Creating a column for LDblock that can be easily converted to gr 
  # Format: chr1:1234-1456
  mutate(grb.coords = str_c(chr,":",start,"-",end)) %>%
  dplyr::select(targets, grb.coords) %>%
  dplyr::rename(GRB_GENES=targets)


# loci.gr is loci.hg19.gr in Anjas script
loci.gr <- loci_in_GRBs.gr
byGRB_list <- vector("list", length(grb_subset))
names(byGRB_list) <- str_replace_all(as.character(grb_subset), "[:-]", "_")
loaded_chr <- "" # for second and third element
for (i in 1:length(grb_subset)) {

  # make each element a list of length 3
  byGRB_list[[i]] <- vector("list", 4)
  print(paste0("Processing ", names(byGRB_list)[[i]], " GRB [", i, "/", length(grb_subset), "]"))
  
  # fill in the first element - df of gene classification for that GRB
  byGRB_list[[i]][[1]] <- per_GRB_genes_classification(grb_subset[i], loci_in_GRBs.gr, genes.hg19.ensembl = genes.hg19.ensembl)
  
  if(!identical(loaded_chr, as.character(seqnames(grb_subset[i])))) {
    
    loaded_chr <- as.character(seqnames(grb_subset[i]))
    #rho_gr <- load_per_chr_rho_file(loaded_chr)
    empP_gr <- load_per_chr_empP_file(loaded_chr, emp_dir = file.path(data_dir, "processed","empP_files_denovo"))
    
  }
  
  # fill in the second element - df of spearman's rho
  #byGRB_list[[i]][[2]] <- per_GRB_correlation_matrix(rho_gr, grb_subset[i], byGRB_list[[i]][[1]]$ens)
  
  # fill in the fourth element with empP of seeing enhancer expression and high promoter expression for the same samples
  
  byGRB_list[[i]][[3]] <- per_GRB_empP_matrix(empP_gr, grb_subset[i], byGRB_list[[i]][[1]]$ens)
    
  
  # fill in the fourth element with information on enhancers in this GRB
  # byGRB_list[[i]][[4]] <- 
  #    per_GRB_enhancers_classification(fantom5.enhancers.tpm.gr, 
  #                                     grb_subset[i],
  #                                     byGRB_list[[i]][[1]]$ens,
  #                                     loci.gr, 
  #                                     TSSs.perTranscript.gr)
}
# Last
saveRDS(byGRB_list, file = file.path(data_dir,"outputs","byGRB_list_plink_updated.rds"))
# saveRDS(loci_in_GRBs.gr, file.path("loci_in_GRBs.gr.rds"))
# Save as outputs table
write.csv(loci_in_GRBs.df, file = file.path(data_dir,"outputs", "loci_in_GRBs.df.csv"))
# Save the rds object to processed
saveRDS(loci_in_GRBs.df, file = file.path(data_dir, "processed","loci_in_GRBs.df.rds"))
# Making GRanges for the same dataframe, now separated
loci_in_GRBs.gr <- grangesFromDataframe(loci_in_GRBs.df, LDblock)
