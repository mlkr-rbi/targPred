# Description: Using empirical approach and geneHancer data to select
# enhancers that are long or short range and extracting sequences

library(GenomicRanges)
library(plyranges)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(rtracklayer)

# sequences for hg19
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)


project_dir <- "/home/katarina/MyProjects/targPred"
inputs_dir <- "data/inputs"
processed_dir <- "data/processed"
outputs_dir <- "/data/outputs"

source(file.path(project_dir, "src/scripts/utils.R"))

# Loading all geneHancer scores data
geneHancer <- read.table(file.path(project_dir, inputs_dir, "geneHancer/GeneHancer_AnnotSV_gene_association_scores_v5.16.txt"),
  header = TRUE)

# Anjas granges object that contains ensembl and hgnc ids
hg19genes <- readRDS(file.path(project_dir, processed_dir, "genes.hg19.ensembl.rds"))

# gh that are from FANTOM5 to reduce number of enh
gh_tissues <- read.table(file.path(project_dir, inputs_dir, "geneHancer/GeneHancer_AnnotSV_tissues_v5.16.txt"),
  header = TRUE, sep = "\t") %>%
  filter(source == "FANTOM5")

# Reducing number of gh
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

data_ <- as.data.frame(gh_empP_ov) %>%
  mutate(gh_pos = paste0(seqnames,":", start,"-", end)) %>%
  select(-seqnames, -start, - end, -strand) %>%
  group_by(gh_pos) %>%
    summarise(min_EP_dist = min(EP_dist),
    max_EP_dist = max(EP_dist),
    enh_count = n(),
    median_EP_dist = median(EP_dist),
    gh_score_mean = mean(combined_score),
    is_elite = median(is_elite),
    eff_size_mean = mean(eff_size),
    r_mean = mean(spearman_r),
    E_active = median(E_active))


data_ <- data_ %>%
  filter(eff_size_mean < 9 ) %>%
  filter(median_EP_dist < 2*10^6)

ggplot(filter(data_, enh_count >=2), aes(x=log10(median_EP_dist), y =gh_score_mean, color = E_active)) +
  geom_point()

print("Checkpoint 2")

x <- data_ %>%
  select(-gh_pos) %>%
  as.data.frame() %>%
  mutate(
    min_EP_dist = scale(min_EP_dist),
    max_EP_dist = scale(max_EP_dist),
    enh_count = scale(enh_count),
    median_EP_dist = scale(median_EP_dist),
    gh_score_mean = scale(gh_score_mean),
    E_active = scale( E_active)) %>% 
  # select(-min_EP_dist) %>%
  # select(-max_EP_dist) %>%
  select(-median_EP_dist) %>%
  select(everything())


# x <- scale(as.numeric(x$min_EP_dist))

cl_num <- 2

k <- kmeans(x=x, centers=cl_num, iter.max = 100, nstart = 100, algorithm = c("Hartigan-Wong"), trace = FALSE)

data_$cluster <- as.factor(k$cluster)

ggplot(data_, aes(x=median_EP_dist, y = enh_count, color = cluster)) +
  geom_point() +
  scale_color_manual(values = brewer.pal(cl_num, "Set1"))


ggplot(data_, aes(x = median_EP_dist, fill = cluster)) +
  geom_density(alpha = 0.5) 

lr_cluster <- ifelse(median(filter(data_, cluster ==1)$median_EP_dist) > median(filter(data_, cluster ==2)$median_EP_dist), 1,2 )
sr_cluster <- ifelse(lr_cluster == 2, 1, 2)

# cutoff_long <- filter(data_, cluster == 1) %>% pull(median_EP_dist) %>% median()
cutoff_long <- filter(data_, cluster == lr_cluster) %>% pull(median_EP_dist) %>% quantile(probs = 0.7)

# cutoff_short <- filter(data_, cluster == 2) %>% pull(median_EP_dist) %>% median()
cutoff_short <- filter(data_, cluster == sr_cluster) %>% pull(median_EP_dist) %>% quantile(probs = 0.15)

ggplot(data_, aes(x = median_EP_dist, fill = cluster)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = cutoff_long, linetype="dashed", color = "purple", size=0.9) +
  geom_vline(xintercept = cutoff_short, linetype="dashed", color = "purple", size=0.9) +
  xlim(0,10**6) +
  scale_color_manual(values = brewer.pal(cl_num, "Set1"))

nucleotide_transformer_gr <- data_ %>%
    mutate(class = case_when(median_EP_dist > cutoff_long & cluster == lr_cluster ~ "long-range",
        median_EP_dist < cutoff_short & cluster == sr_cluster ~ "short-range",
        TRUE ~ "not-defined"))

# write.table(nucleotide_transformer_gr, file.path(project_dir, outputs_dir, "nucleotide_transformer/enhancers_annotation_all_gr.txt"), quote=FALSE,row.names=FALSE, sep = "\t")


nucleotide_transformer_gr <- grangesFromDataframe(nucleotide_transformer_gr, colName = gh_pos)

# Filtering out intermediates 
nucleotide_transformer_gr <- nucleotide_transformer_gr[nucleotide_transformer_gr$class != "not-defined"]
# write.table(as.data.frame(mcols(nucleotide_transformer_gr)), file.path(project_dir, outputs_dir, "nucleotide_transformer/enhancers_annotation_short_long_gr.txt"), quote=FALSE,row.names=FALSE, sep = "\t")

# gathering sequences
# write.table(as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg19, nucleotide_transformer_gr)), file.path(project_dir, outputs_dir, "nucleotide_transformer/enh_sequences_short_long.hg19.txt"))

# import.chain does not work and it seems to be the problem with the version so using Genome browser tool
# Download tool with
## wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
# Change mod
## chmod +x liftOver
# Downloading chain file for liftover https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
# Prepare bed file for liftOver
# makeBed(nucleotide_transformer_gr, modify_start = -1, filename= file.path(project_dir, outputs_dir, "nucleotide_transformer/enh_sequences_short_long.hg19.bed"))

# Doing liftOver
## ./liftOver enh_sequences_short_long.hg19.bed hg19ToHg38.over.chain.gz enh_sequences_short_long.hg38.bed unmapped.bed
# unmapped.bed is empty

# nucleotide_transformer_hg38_gr <- rtracklayer::import.bed(file.path(project_dir, outputs_dir, "nucleotide_transformer/enh_sequences_short_long.hg38.bed"))
# mcols(nucleotide_transformer_hg38_gr) <- mcols(nucleotide_transformer_gr)[,-1] # removing hg19 gh_positions
# write.table(as.data.frame(mcols(nucleotide_transformer_hg38_gr)), file.path(project_dir, outputs_dir, "nucleotide_transformer/enhancers_annotation_short_long_gr.hg38.txt"), quote=FALSE,row.names=FALSE, sep = "\t")
# write.table(as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg38, nucleotide_transformer_hg38_gr)), file.path(project_dir, outputs_dir, "nucleotide_transformer/enh_sequences_short_long.hg38.txt"))
