# Description: Preparing an object that has combined information from our empP analysis and geneHancer data

library(GenomicRanges)
library(plyranges)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(mclust, help, pos = 2, lib.loc = NULL)

# sequences for hg19
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(BSgenome.Hsapiens.UCSC.hg38)


project_dir <- "/home/katarina/MyProjects/targPred"
inputs_dir <- "data/inputs"
processed_dir <- "data/processed"
outputs_dir <- "/data/outputs"

source(file.path(project_dir, "src/scripts/utils.R"))

# Enhancer data was requested through a website and sent by email. Everything is uploaded on Orthus in data/inputs/geneHnacer

# Loading all geneHancer scores data
geneHancer <- read.table(file.path(project_dir, inputs_dir, "geneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.16.txt"),
  header = TRUE)

# Anjas granges object that contains ensembl and hgnc ids
hg19genes <- readRDS(file.path(project_dir, processed_dir, "genes.hg19.ensembl.rds"))

# gh that are from FANTOM5 to reduce number of enh
gh_tissues <- read.table(file.path(project_dir, inputs_dir, "geneHancer/GeneHancer_AnnotSV_tissues_v5.16.txt"),
  header = TRUE, sep = "\t") %>%
  filter(source == "FANTOM5")

# Reducing number of gh by selecting enhancers sourced from F5
geneHancer <- geneHancer %>%
  filter( is.element(GHid, gh_tissues$GHid)) %>%
  filter(is.element(symbol, hg19genes$hgnc_symbol))

# Adding ensembl ids
temp <- data.frame(hgnc_symbol = hg19genes$hgnc_symbol, ensembl = names(hg19genes))

geneHancer <-
  dplyr::left_join(geneHancer,temp, by = c("symbol" = "hgnc_symbol"))

# Loading genomic location of enhancers
gh_bed <- read.table(file.path(project_dir, inputs_dir, "geneHancer/GeneHancer_hg19_v5.16.bed"))

gh_bed <- filter(gh_bed, is.element(V4, geneHancer$GHid))

gh_bed <- left_join(gh_bed, geneHancer, by = c("V4" = "GHid"))

gh_gr <- GRanges(seqnames = gh_bed$V1,
  ranges = IRanges(start = gh_bed$V2,
  end = gh_bed$V3))

gh_gr$gh_names <- gh_bed$V4
values(gh_gr) <- select(gh_bed, symbol,ensembl, combined_score, is_elite)

# Genome wide empP data with added metrics
chr_vec <- paste0("chr", c(1:22)) %>%
  append("chrX")

genome_EP_distance <- data.frame()

print("Checkpoint 1")
for(chr in chr_vec){
  
      empP_gr <- read.table(file = file.path(project_dir, processed_dir, "empP_et_all_files_denovo", paste0("empP_CAGEen_denovo_hg19Ens_6.4Mb_", chr, ".txt")),
      header = TRUE)
  
  
  data_ <- as.data.frame(empP_gr) %>%
  as.data.frame() %>%
  dplyr::filter(empP < 0.001)

  genome_EP_distance <- rbind(genome_EP_distance, data_)
}


# dropping chromosome column since it causes errors
genome_EP_distance <- genome_EP_distance %>%
  select (-chr)
genome_EP_distance_gr <- grangesFromDataframe(genome_EP_distance, colName = enh_pos)

gh_empP_ov <- join_overlap_intersect(gh_gr, genome_EP_distance_gr)
gene_ids <- gh_empP_ov$gene
gene_ids <- str_extract(gene_ids, "^[^:]+")
gh_empP_ov <- plyranges::filter(gh_empP_ov, ensembl == gene_ids)

# Saving gh_empP_ov as rds
saveRDS(gh_empP_ov, file.path(project_dir, processed_dir, "geneHancer_fantom5_overlap.gr.rds"))
