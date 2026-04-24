### Description ####
# Combining GRBs from bed file and target genes to make a GRanges 


library(tidyverse)
library(plyranges)
library(GenomicRanges)
library(dplyr)
library(rtracklayer)

path_file <- "src/core/paths.txt"
path_vector <- readr::read_lines(path_file)
directory_type <- ""
data_dir <- ""

for (idx in 1:length(path_vector)) {
    line_vector <- strsplit(path_vector[idx],"")
    if (line_vector[[1]][1] == "$") {
        directory_type <- paste(line_vector[[1]][-1], collapse = "")
        next
    } else {
        if (directory_type == "data") data_dir <- path_vector[idx]
    }
}

if (data_dir == "") {
    print("ERROR: data_dir not defined in paths.txt!")
    quit(status=1)
}

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

### Load data ####
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


# Building target gene GRanges
## Removing target genes from the patched sequences
target_genes <- filter(target_genes, !str_detect(ucscCoordinates.hg19, "PATCH"))

target_genes_gr <- grangesFromDataframe(target_genes,
                                        colName =ucscCoordinates.hg19)

# Adding target genes to GRB and 
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

# Add number of target genes to run_targPred log file
num_target_genes <- sum(cne_density_grb$geneName != "")
log_file <- file.path(data_dir, "processed/logs", "targPred_run_log.txt")
write(paste("Number of target genes in CNE density GRBs:", num_target_genes), file = log_file, append = TRUE)
# Add number of GRBs 
num_grbs <- length(cne_density_grb)
write(paste("Number of GRBs with target genes:", num_grbs), file = log_file, append = TRUE)
# Add number of GRBs with multiple target genes
num_grbs_multiple_targets <- sum(str_count(cne_density_grb$geneName, ",") >= 1)
write(paste("Number of GRBs with multiple target genes:", num_grbs_multiple_targets), file = log_file, append = TRUE)



# Save as granges object
saveRDS(cne_density_grb, file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))
# Save as table for use in scripts that need a table
cne_density_grb_df <- as.data.frame(cne_density_grb) %>%
  select(-strand)

write.table(cne_density_grb_df, file.path(data_dir, "processed", "GRBs_with_targets.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)  

