#### Libraries #####
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
# library(rols)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)


source(file.path("src/core", "utils.R"))

path_file <- "src/core/paths.txt"
paths <- set_paths(path_file)
core_dir <- paths$core_dir
work_dir <- paths$work_dir
data_dir <- paths$data_dir
empP_dir <- paths$empP_dir
gwas_year <- paths$gwas_year
LDblock_path <- paths$LDblock_path
byGRBlist_name <- paths$byGRBlist_name



setwd(work_dir)


# Sent by Anja by mail (from fantom5 data)
TSSs.perTranscript.gr <- readRDS(file.path(data_dir, "inputs","TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(data_dir, "inputs","fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(data_dir, "inputs","promoter_cage_tpm_df_fixed_pooled.rds"))


# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(data_dir, "inputs","genes.hg19.ensembl.rds"))

# Load GRBs with target genes
grbs <- readRDS(file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))

# Load loci in GRBs as dataframe and convert to GRanges for later use
loci_in_GRBs.df <- read.csv(file.path(data_dir, "processed", "loci_in_GRBs.df.csv"), stringsAsFactors = FALSE) 
# FIlter rows that have NA in MAPPED_TRAIT column
loci_in_GRBs.df <- loci_in_GRBs.df %>%
  filter(!is.na(MAPPED_TRAIT))


# Making GRanges for the same dataframe, now separated
loci_in_GRBs.gr <- grangesFromDataframe(loci_in_GRBs.df, LDblock) 

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



byGRB_list <- vector("list", length(grb_subset))
names(byGRB_list) <- str_replace_all(as.character(grb_subset), "[:-]", "_")
loaded_chr <- "" # for second and third element

# empP_gr <- load_per_chr_empP_file("chr1", emp_dir = file.path(data_dir, "processed",empP_dir))

for (i in 1:length(grb_subset)) {

    # make each element a list of length 3
    byGRB_list[[i]] <- vector("list", 4)
    print(paste0("Processing ", names(byGRB_list)[[i]], " GRB [", i, "/", length(grb_subset), "]"))

    # Paste number of GRBs processed to log file
    cat(paste0("Processing ", names(byGRB_list)[[i]], " GRB [", i, "/", length(grb_subset), "]\n"), file = file.path(data_dir, "processed/logs", "byGRB_list_processing.txt"), append = TRUE)

    # fill in the first element - df of gene classification for that GRB
    byGRB_list[[i]][[1]] <- per_GRB_genes_classification(grb_subset[i], loci_in_GRBs.gr, genes.hg19.ensembl = genes.hg19.ensembl)
    
    if(!identical(loaded_chr, as.character(seqnames(grb_subset[i])))) {

        print(file.path(data_dir, "processed",empP_dir))
      
        loaded_chr <- as.character(seqnames(grb_subset[i]))
        empP_gr <- load_per_chr_empP_file(loaded_chr, emp_dir = file.path(data_dir, "processed",empP_dir))
      
    }

    # Paste size of the gene classification df to log file
    cat(paste0("Size of gene classification df for ", names(byGRB_list)[[i]], " GRB: ", nrow(byGRB_list[[i]][[1]]), " genes\n"), file = file.path(data_dir, "processed/logs", "byGRB_list_processing.txt"), append = TRUE)
    
    # fill in the second element - df of spearman's rho
    #byGRB_list[[i]][[2]] <- per_GRB_correlation_matrix(rho_gr, grb_subset[i], byGRB_list[[i]][[1]]$ens)
    
    # fill in the fourth element with empP of seeing enhancer expression and high promoter expression for the same samples
    
    byGRB_list[[i]][[3]] <- per_GRB_empP_matrix(empP_gr, grb_subset[i], byGRB_list[[i]][[1]]$ens, loci.gr = loci_in_GRBs.gr)
    
    # paste number of genes with empP values for that GRB to log file
    cat(paste0("Number of genes with empP values for ", names(byGRB_list)[[i]], " GRB: ", sum(!is.na(byGRB_list[[i]][[3]]$empP)), " out of ", nrow(byGRB_list[[i]][[3]]), " genes\n"), file = file.path(data_dir, "processed/logs", "byGRB_list_processing.txt"), append = TRUE)
    
    # fill in the fourth element with information on enhancers in this GRB
    # byGRB_list[[i]][[4]] <- 
    #    per_GRB_enhancers_classification(fantom5.enhancers.tpm.gr, 
    #                                     grb_subset[i],
    #                                     byGRB_list[[i]][[1]]$ens,
    #                                     loci.gr, 
    #                                     TSSs.perTranscript.gr)
}
### Saving data ####
# write.csv(gwas, file = file.path(data_dir,"outputs", "gwas_associations_LD_processed.csv"))
saveRDS(byGRB_list, file = file.path(data_dir,"processed", byGRBlist_name))
# Save as outputs table
# write.csv(loci_in_GRBs.df, file = file.path(data_dir,"outputs", "loci_in_GRBs.df.csv"))
# Save the rds object to processed
# saveRDS(loci_in_GRBs.df, file = file.path(data_dir, "processed","loci_in_GRBs.df.rds"))
# Making GRanges for the same dataframe, now separated
# loci_in_GRBs.gr <- grangesFromDataframe(loci_in_GRBs.df, LDblock)
# saveRDS(loci_in_GRBs.gr, file = file.path(data_dir, "processed","loci_in_GRBs.gr.rds"))
library(biomaRt)
library(GenomicRanges)

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl",
  host = "https://grch37.ensembl.org"
)

attrs <- c(
  "chromosome_name",
  "transcript_start",
  "transcript_end",
  "strand",
  "ensembl_gene_id",
  "ensembl_transcript_id",
  "hgnc_symbol"
)

bm <- getBM(attributes = attrs, mart = mart)
bm <- bm[bm$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
bm$seqnames <- paste0("chr", bm$chromosome_name)

gr <- GRanges(
  seqnames = bm$seqnames,
  ranges = IRanges(start = as.integer(bm$transcript_start),
                   end = as.integer(bm$transcript_end)),
  strand = ifelse(bm$strand == 1, "+", "-"),
  ensembl_gene_id = bm$ensembl_gene_id,
  transcript_id = bm$ensembl_transcript_id,
  hgnc_symbol = bm$hgnc_symbol
)

saveRDS(gr, "data/processed/gene_model_hg19.gr.rds")
