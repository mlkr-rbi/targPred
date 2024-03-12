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

#### Libraries #####
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(rols)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

path_file <- "paths.txt"
path_vector <- readr::read_lines("paths.txt")
directory_type <- ""
core_dir <- ""
work_dir <- ""
data_dir <- ""
empP_dir <- ""
gwas_year <- ""
LDblock_path <- ""
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
        if (directory_type == "empP_dir") empP_dir <- path_vector[idx]
        if (directory_type == "LDblock_path") LDblock_path <- path_vector[idx]
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
} else if (empP_dir == "") {
    print("ERROR: empP_dir not defined in paths.txt!")
    quit(status=1)
} else if (gwas_year == "") {
    print("ERROR: gwas_year not defined in paths.txt!")
    quit(status=1)
}
if (LDblock_path == "") {
    print("ERROR: full path to bed file with LD blocks not defined in paths.txt!")
    quit(status=1)
}
if (byGRBlist_name == "") {
    print("byGRB_list name is default byGRB_list.rds")
    byGRBlist_name <- "byGRB_list.rds"
}

source(file.path(core_dir, "utils.R"))

setwd(work_dir)

# create GWAS associations here instead of a separate Python script

gwas_associations <- readRDS(file.path(data_dir, "inputs", paste0("gwas_catalog_",gwas_year),"gwas_associations.rds"))
gwas_snps_LDblocks = read.csv(file.path(LDblock_path),sep="\t",header=FALSE)
colnames(gwas_snps_LDblocks) <- c("chromosome","min_loc","max_loc","SNPS")
gwas_snps_LDblocks <-
  gwas_snps_LDblocks %>%
  # Creating a column for LDblock that can be easily converted to gr 
  # Format: chr1:1234-1456
  mutate(LDblock = str_c(chromosome,":",min_loc,"-",max_loc)) %>%
  dplyr::select(SNPS, LDblock)

gwas <- left_join(gwas_snps_LDblocks, gwas_associations, by="SNPS")
saveRDS(gwas , file.path(data_dir, "outputs", "gwas_associations_LDblock.rds"))

grbs <- readRDS(file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))
 

# Target genes for CNE density GRBs 
target_genes <- 
    read.csv(file.path(data_dir, "inputs","hg19.canFam3.tsv"), sep = "\t", header = TRUE)

# Downloaded efo terms for trait separating
efo_terms <- readRDS(file.path(data_dir, "inputs","EFO_terms.rds"))

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
    data_r <- makeGRangesFromDataFrame(data_r,
                                       seqnames.field ="chr",
                                       start.field = "start",
                                       end.field = "end",
                                       keep.extra.columns = keep_metadata)
    return(data_r)
}

 
### Modifying genes from gwas data #####
# Gwas catalog stores the putative target gene of an SNP in two columns: (from https://www.ebi.ac.uk/gwas/docs/fileheaders ) REPORTED GENE(S): Gene(s) reported by author MAPPED GENE(S): Gene(s) mapped to the strongest SNP. If the SNP is located within a gene, that gene is listed. If the SNP is located within multiple genes, these genes are listed separated by commas. If the SNP is intergenic, the upstream and downstream genes are listed, separated by a hyphen. - Combine the two columns in a single column GWAS_genes

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
    mutate(temp1 = map_chr(str_split(temp, ", |;| - |-"), ~str_c(unique(.x),collapse = ", "))) %>%
    mutate(GWAS_GENES = map_chr(str_split(temp1, ", "), ~str_c(.x[!.x == "NA"],collapse = ", "))) %>%
    dplyr::select(-c(temp, temp1, REPORTED.GENE.S., MAPPED_GENE))


### Separating Traits ####
# Each investigated SNP has an associated mapped trait, or multiple. Since our goal is to explore the relationship between a SNP and a single trait row separation based on trait is required. Trait names are stored in MAPPED_TRAIT and the URI that contain unique disease/trait ontology id in MAPPED_TRAIT_URI. Traits are comma separated but some trait names contain commas. To avoid error prone case by case exceptions, trait separation will be based on the URI and trait names will be cross mapped from the ontology.

## Trait ontology contains all GWAS trait EFO IDs extracted from MAPPED_TRAIT_URI
trait_ontology <- dplyr::select(gwas, MAPPED_TRAIT_URI) %>% unique()
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

### Relabeling Traits ####
# Using trait_ontology dataframe to do this once for each term
trait_ontology$label <- NA

for(id in 1:length(trait_ontology$EFO_ID)){
  
    term_id <- str_replace(trait_ontology$EFO_ID[id], "_",":")
    
    # efo_terms is an ontology object which needs the termLabel method to access label
    tryCatch(trait_ontology$label[id] <- rols::termLabel(efo_terms[term_id]), error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

gwas <-
    gwas %>%
    # Multiple URIs that correspond to one trait are sep with ","
    separate_rows(MAPPED_TRAIT_URI, sep = ",") %>%
    mutate(MAPPED_TRAIT_URI = stringr::str_trim(MAPPED_TRAIT_URI, side = "both")) %>%
    # Splitting URI slashes since the unique ID is stored in the end
    mutate(MAPPED_TRAIT_ID = str_split(MAPPED_TRAIT_URI, pattern = "/", simplify =TRUE)[,5]) %>%
    left_join(dplyr::select(trait_ontology, MAPPED_TRAIT_URI, label), by = ("MAPPED_TRAIT_URI")) %>%
    filter(MAPPED_TRAIT_URI != "")

# Removing the obsolete_ addition from some trait labels
gwas <-
    gwas %>%
    mutate(MAPPED_TRAIT = label) %>%
    dplyr::select(-label) %>%
    mutate(MAPPED_TRAIT = ifelse(str_detect(MAPPED_TRAIT,"obsolete_"),
                                 str_remove(MAPPED_TRAIT, "obsolete_"), 
                                 MAPPED_TRAIT)) %>%
    mutate(MAPPED_TRAIT = ifelse(str_detect(MAPPED_TRAIT,"obsolete "),
                                 str_remove(MAPPED_TRAIT, "obsolete "), 
                                 MAPPED_TRAIT)) %>%
    # Removing NA MAPPED_TRAIT
    filter(!is.na(MAPPED_TRAIT))
  
write.csv(gwas, file = file.path(data_dir,"outputs", "gwas_with_traits.csv"))

## gwas_gr stores SNP LD blocks in GenomicRanges object
gwas.gr <- grangesFromDataframe(gwas, LDblock)

## Loci in GRBs ####

# Overlaping gwas loci with GRBs
loci_in_GRBs.gr <- subsetByOverlaps(gwas.gr, grbs)
ovlp <- findOverlaps(loci_in_GRBs.gr, grbs)

# Adding GRB coordinates to each locus as character
loci_in_GRBs.gr$grb.coords[queryHits(ovlp)] <-
    as.character(grbs[subjectHits(ovlp)])

# Adding GRB target genes 
loci_in_GRBs.gr$GRB_GENES[queryHits(ovlp)] <- 
    grbs[subjectHits(ovlp)]$geneName

# Saving loci in GRB as bed file for visualizations in Genome browser
export.bed(loci_in_GRBs.gr, con = file.path(data_dir, "outputs", "loci_in_GRBs.bed"))

loci_in_GRBs.df <- mcols(loci_in_GRBs.gr) %>% as.data.frame()

# SNPs in GRBs granges

# Coercing to data frame because 
# it makes a GPos object by default which makes bugs afterwards
snps_in_GRBs.gr <- 
  snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh37, 
           unique(loci_in_GRBs.df$SNPS),
           ifnotfound = "drop") %>%
  as.data.frame() %>%
  dplyr::select(seqnames, pos, SNPS = "RefSNP_id")
# # Generating Granges 
snps_in_GRBs.gr <- 
  makeGRangesFromDataFrame(snps_in_GRBs.gr,
                           seqnames.field = "seqnames",
                           start.field = "pos",
                           end.field ="pos",
                           keep.extra.columns = TRUE)


# Concetinating  to the final object
## Changing type of seqnames from 1,2,3,... to chr1, chr2, chr3,...
 seqlevelsStyle(snps_in_GRBs.gr) <- "UCSC"
# 
# # ?Maybe combine with loci_in_GRBs.df for extra metadata
 saveRDS(snps_in_GRBs.gr, file = file.path(data_dir, "processed","snps_in_GRBs.gr.rds"))

####### byGRB_list ##########

# Sent by Anja by mail (from fantom5 data)
TSSs.perTranscript.gr <- readRDS(file.path(data_dir, "processed","TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(data_dir, "processed","fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(data_dir, "processed","promoter_cage_tpm_df_fixed_pooled.rds"))


# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(data_dir, "processed","genes.hg19.ensembl.rds"))


# makes sure to work on the GRB level, even when one GRB covers multiple loci
grb_subset <- subsetByOverlaps(grbs, loci_in_GRBs.gr)

# add grb targets to loci in grbs
grb_subset.df <- data.frame(chr=seqnames(grb_subset),
                            start=start(grb_subset),
                            end=end(grb_subset),
                            targets=grb_subset$geneName,
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
        empP_gr <- load_per_chr_empP_file(loaded_chr, emp_dir = file.path(data_dir, "processed",empP_dir))
      
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
### Saving data ####
write.csv(gwas, file = file.path(data_dir,"outputs", "gwas_associations_LD_processed.csv"))
saveRDS(byGRB_list, file = file.path(data_dir,"outputs",byGRBlist_name))
# Save as outputs table
write.csv(loci_in_GRBs.df, file = file.path(data_dir,"outputs", "loci_in_GRBs.df.csv"))
# Save the rds object to processed
saveRDS(loci_in_GRBs.df, file = file.path(data_dir, "processed","loci_in_GRBs.df.rds"))
# Making GRanges for the same dataframe, now separated
loci_in_GRBs.gr <- grangesFromDataframe(loci_in_GRBs.df, LDblock)
saveRDS(loci_in_GRBs.gr, file = file.path(data_dir, "processed","loci_in_GRBs.gr.rds"))
