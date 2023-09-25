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


library(GenomicRanges)
library(rtracklayer)
library(purrr)
library(curl)
library(RColorBrewer)
library(org.Hs.eg.db)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(tidyverse)
library(httr)
library(jsonlite)
library(xml2)

source(file.path("utils.R"))

aux_dir <- "../aux_files"

# Output from LD_Block_Calculation.Rmd 
 gwas <- read.csv("../gwas_associations_LDblock.csv") %>%
   select(CHR= "CHR_ID", SNPS, LDblock, STUDY.ACCESSION, REPORTED.GENE.S., 
          MAPPED_GENE, P.VALUE, MAPPED_TRAIT, MAPPED_TRAIT_URI) %>%
   filter(!is.na(CHR))
# Main GRB dataset
cne_density_grb <- 
  import.bed(file.path(aux_dir,"hg19.canFam3.98.50.processed.threshold.grbs.bed")) %>%
  makeGRangesFromDataFrame() %>%
  # Reducing since some GRBs overlap
  # Look at CNE_density_overlapping_GRBs.txt for examples
  GenomicRanges::reduce()
start(cne_density_grb) = start(cne_density_grb) -1
# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(aux_dir,"hg19.canFam3.tsv"), sep = "\t", header = TRUE)



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
  # Creating a temporary column that contains all genes comma sep
  mutate(temp = case_when(is.element(REPORTED.GENE.S., non_gene_names) ~ MAPPED_GENE,
                                is.element(MAPPED_GENE, non_gene_names) ~ REPORTED.GENE.S.,
                                TRUE ~ str_c(REPORTED.GENE.S., MAPPED_GENE, sep = ", "))) %>%
  # Removing repeated genes (eg. SKIV2L, SKIV2L, SKIV2L, SKIV2L, SKIV2L, SKIV2L)
  mutate(GWAS_GENES = map_chr(str_split(temp, ", "), ~str_c(unique(.x),collapse = ", "))) %>%
  dplyr::select(-c(temp, REPORTED.GENE.S., MAPPED_GENE))
## gwas_gr stores SNP LD blocks in GenomicRanges object
gwas_gr <- grangesFromDataframe(gwas, LDblock)


# Overlapping all SNps to GRBs(After modification there should not be duplicated query hits)
ov_gwas_cne <-
  findOverlaps(query = gwas_gr, subject = cne_density_grb)
# Around 500 SNPs that overlap two GRBs
duplicated(queryHits(ov_gwas_cne)) %>% table()
# Storing overlap results in a dataframe and selecting only hits that correspond to SNP LD block that overlaps two GRBs
snps_two_grbs_index <- 
  as.data.frame(ov_gwas_cne) %>%
  filter(is.element(queryHits, queryHits(ov_gwas_cne)[duplicated(queryHits(ov_gwas_cne))]))
# Storing granges of SNPs that span two GRBs for visualization
snps_two_grbs <- gwas_gr[snps_two_grbs_index$queryHits,]
# Commenting out this so that every rerun does not overwrite
export.bed(snps_two_grbs, "SNPs_two_GRBs.bed")
rm(snps_two_grbs, snps_two_grbs_index)


# Building target gene genomic ranges
target_genes_gr <- grangesFromDataframe(target_genes,
                                        colName =ucscCoordinates.hg19)
# Adding target genes as metadata 
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


loci_in_GRBs.gr <- subsetByOverlaps(gwas_gr, cne_density_grb)
ovlp <- findOverlaps(loci_in_GRBs.gr, cne_density_grb)
# Adding GRB coordinates to each locus as character
loci_in_GRBs.gr$grb.coords[queryHits(ovlp)] <-
  as.character(cne_density_grb[subjectHits(ovlp)])

# Adding GRB target genes 
loci_in_GRBs.gr$GRB_GENES[queryHits(ovlp)] <- 
  cne_density_grb[subjectHits(ovlp)]$geneName

# Loci in GRBs

temp <-
  data.frame(seqnames=seqnames(loci_in_GRBs.gr),
  starts=start(loci_in_GRBs.gr)-1,
  ends=end(loci_in_GRBs.gr),
  names=c(rep(".", length(loci_in_GRBs.gr))),
  scores=c(rep(".", length(loci_in_GRBs.gr))),
  strands=strand(loci_in_GRBs.gr))
write.table(temp, file="loci_in_GRBs.bed", quote=F, sep="\t", row.names=F, col.names=F)
export.bed(loci_in_GRBs.gr, con = "loci_in_GRBs.bed")

# Separating Traits

loci_in_GRBs.df <- mcols(loci_in_GRBs.gr) %>%
  as.data.frame()
trait_ontology <- dplyr::select(loci_in_GRBs.df, MAPPED_TRAIT_URI) %>%
  unique()
rownames(trait_ontology) <- NULL
trait_ontology <- 
  trait_ontology %>%
  # Separating multiple traits into rows
  tidyr::separate_rows(MAPPED_TRAIT_URI, sep = ",") %>%
  unique() %>%
  # Splitting the URI eg. http://www.ebi.ac.uk/efo/EFO_0007978 to select only the ID EFO_0007978
  mutate(EFO_ID = str_split(MAPPED_TRAIT_URI, pattern = "/", simplify = TRUE)[,5]) %>%
  filter(!is.na(MAPPED_TRAIT_URI))


## Downloading EFO ontology

# Downloaded efo terms
efo_terms <- readRDS(file.path(aux_dir,"EFO_terms.rds"))
# Initializing empty df
efo_mappings <- data.frame(id1=NA, id2=NA, label=NA) 
# Looping through elements of the list
for(i in 1:length(names(efo_terms@x))){
  id <- names(efo_terms@x)[i]
  if(i%%100 == 0) {message(i)}
    if(is.null(unlist(efo_terms[[id]]@annotation$database_cross_reference))){
      # Some ids have no cross refrences
      working_df <-  data.frame(id, NA, efo_terms[[id]]@label)
    }
  
  else{
  working_df <-  data.frame(id, unlist(efo_terms[[id]]@annotation$database_cross_reference), efo_terms[[id]]@label)
  }
  
  # Changing colnames of working_df before binding
  colnames(working_df) <- c("id1","id2","label")
  efo_mappings <- rbind(efo_mappings, working_df)
}
efo_mappings <- efo_mappings %>%
  # Removing a lot of duplicates
  unique() %>%
  # Removing first row that is only NA
  filter(!is.na(id1) & !is.na(id2) & !is.na(label))
write.csv(efo_mappings, file = "EFO_mappings.csv")

## Relabeling Traits

# Loading efo_terms
#efo_terms <- readRDS(file.path(aux_dir,"EFO_terms.rds"))
loci_in_GRBs.df <-
  loci_in_GRBs.df %>%
  # Multiple URIs that correspond to one trait are sep with ","
  separate_rows(MAPPED_TRAIT_URI, sep = ",") %>%
  # Splitting URI slashes since the unique ID is stored in the end
  mutate(MAPPED_TRAIT_ID = str_split(MAPPED_TRAIT_URI, pattern = "/", simplify =TRUE)[,5]) %>%
  # Replacing underscore "_" with colon ":" since it is standard for ontology ids
  mutate(MAPPED_TRAIT_ID = str_replace_all(MAPPED_TRAIT_ID, "_",":")) 
for(id in 1:dim(loci_in_GRBs.df)[1]){
  
  term_id <- loci_in_GRBs.df$MAPPED_TRAIT_ID[id]
  # efo_terms is an ontology object which needs the termLabel method to access label
  tryCatch(loci_in_GRBs.df$MAPPED_TRAIT[id] <- rols::termLabel(efo_terms[term_id]),
           error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}
# Removing the obsolete_ addition from some trait labels
loci_in_GRBs.df <- 
  loci_in_GRBs.df %>%
  mutate(MAPPED_TRAIT = ifelse(str_detect(MAPPED_TRAIT,"obsolete_"),
                               str_remove(MAPPED_TRAIT, "obsolete_"), 
                               MAPPED_TRAIT)) %>%
  # Removing NA MAPPED_TRAIT
  filter(!is.na(MAPPED_TRAIT))
  
# Save as outputs table
write.csv(loci_in_GRBs.df, file = file.path(project_dir, outputs_dir, "loci_in_GRBs.df.csv"))
# Save the rds object to processed
saveRDS(loci_in_GRBs.df, file = file.path("loci_in_GRBs.df.rds"))
# Making GRanges for the same dataframe, now separated
loci_in_GRBs.gr <- grangesFromDataframe(loci_in_GRBs.df, LDblock)


# SNPs in GRBs granges

# Coercing to data frame because 
# it makes a GPos object by default which makes bugs afterwards
snps_in_GRBs.gr <- 
  snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, 
           loci_in_GRBs.df$SNPS %>% unique(),
           ifnotfound = "drop") %>%
  as.data.frame() %>%
  dplyr::select(seqnames, pos, SNPS = "RefSNP_id")
# Generating Granges 
snps_in_GRBs.gr <- 
  makeGRangesFromDataFrame(snps_in_GRBs.gr,
                           seqnames.field = "seqnames",
                           start.field = "pos",
                           end.field ="pos",
                           keep.extra.columns = TRUE)


# Concetinating  to the final object
#ovo sam ja izkomentirao
#snps_in_GRBs.gr <- c(snps_in_GRBs.gr, missing_snps.gr) %>%
#  unique()
# Changing type of seqnames from 1,2,3,... to chr1, chr2, chr3,...
seqlevelsStyle(snps_in_GRBs.gr) <- "UCSC"
# ?Maybe combine with loci_in_GRBs.df for extra metadata
saveRDS(snps_in_GRBs.gr, file = file.path("snps_in_GRBs.gr.rds"))

# byGRB_list
#loci_in_GRBs.df <- readRDS(file.path("loci_in_GRBs.df.rds"))
# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(aux_dir,"hg19.canFam3.tsv"), sep = "\t", header = TRUE)  
# Building target gene genomic ranges
target_genes_gr <- grangesFromDataframe(target_genes, ucscCoordinates.hg19)
# Sent by Anja by mail (from fantom5 data)
TSSs.perTranscript.gr <- readRDS(file.path(aux_dir,"TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(aux_dir,"fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(aux_dir,"promoter_cage_tpm_df_fixed_pooled.rds"))
#grbs <- readRDS(file.path("GRBs_CNE_density_modf_kurtosis.rds"))
grbs <- cne_density_grb
# Adding target genes to GRBs 
grbs$targets <- ""
for (i in 1:length(grbs)) {
  t <- sort(subsetByOverlaps(target_genes_gr, grbs[i]))
  if(length(t)>0) {
    grbs$targets[i] <- paste(t$geneName, collapse=",")
  }
}
 saveRDS(grbs, file.path("GRBs_with_targets.gr.rds"))
# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(aux_dir,"genes.hg19.ensembl.rds"))
gwas_associations <- readRDS(file.path(aux_dir,"gwas_associations.rds"))
# Started building snps in GRB not fiinished!!
snps_in_GRBs.gr <-
  as.data.frame(mcols(loci_in_GRBs.gr)) %>%
  left_join(dplyr::select(gwas_associations, SNPS, CHR_POS))


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

saveRDS(loci_in_GRBs.df, file = file.path("loci_in_GRBs.df.rds"))
saveRDS(grb_subset, "pop_EUR_2022_grb_subset_GeThesisTargets.rds")
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
  byGRB_list[[i]][[1]] <- per_GRB_genes_classification(grb_subset[i], loci_in_GRBs.gr)
  
  if(!identical(loaded_chr, as.character(seqnames(grb_subset[i])))) {
    
    loaded_chr <- as.character(seqnames(grb_subset[i]))
    #rho_gr <- load_per_chr_rho_file(loaded_chr)
    empP_gr <- load_per_chr_empP_file(loaded_chr, emp_dir = file.path(aux_dir,"empP_files"))
    
  }
  
  # fill in the second element - df of spearman's rho
  #byGRB_list[[i]][[2]] <- per_GRB_correlation_matrix(rho_gr, grb_subset[i], byGRB_list[[i]][[1]]$ens)
  
  # fill in the fourth element with empP of seeing enhancer expression and high promoter expression for the same samples
  byGRB_list[[i]][[3]] <- per_GRB_empP_matrix(empP_gr, grb_subset[i], byGRB_list[[i]][[1]]$ens)
    
  # fill in the fourth element with information on enhancers in this GRB
  byGRB_list[[i]][[4]] <- 
     per_GRB_enhancers_classification(fantom5.enhancers.tpm.gr, 
                                      grb_subset[i],
                                      byGRB_list[[i]][[1]]$ens,
                                      loci.gr, 
                                      TSSs.perTranscript.gr)
}
# Last
saveRDS(byGRB_list, file = file.path("byGRB_list.rds"))
saveRDS(loci_in_GRBs.gr, file.path("loci_in_GRBs.gr.rds"))
