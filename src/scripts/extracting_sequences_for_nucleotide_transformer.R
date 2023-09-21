# Description: Using empirical approach and geneHancer data to select
# enhancers that are long or short range and extracting sequences
# Tried out gmm and kmeans, and kmeans gave better results

library(GenomicRanges)
library(plyranges)
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(rtracklayer)
library(GGally, help, pos = 2, lib.loc = NULL)
library(mclust, help, pos = 2, lib.loc = NULL)

# sequences for hg19
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)


project_dir <- "/home/katarina/MyProjects/targPred"
inputs_dir <- "data/inputs"
processed_dir <- "data/processed"
outputs_dir <- "/data/outputs"

source(file.path(project_dir, "src/scripts/utils.R"))

gh_empP_ov <- readRDS(file.path(project_dir, processed_dir, "geneHancer_fantom5_overlap.gr.rds"))

# Summarising by enhancer
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
    E_active = mean(E_active),
    empP_median = median(empP),
    width = unique(width))


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
    width = scale(width),
    median_EP_dist = scale(median_EP_dist),
    gh_score_mean = scale(gh_score_mean),
    E_active = scale( E_active)) %>% 
  # select(-min_EP_dist) %>%
  select(-max_EP_dist) %>%
  # select(-median_EP_dist) %>%
  # select(-cluster) %>%
  select(everything())


timestamp <- Sys.time()
timestamp <- format(timestamp, "%Y%m%d_%H%M_")



BIC <- mclustBIC(x)
plot(BIC)
summary(BIC)

cl_num <- 3

# gmm <- Mclust(x, G = cl_num, modelName = "VEV")
k <- kmeans(x=x, centers=cl_num, iter.max = 100, nstart = 100, algorithm = c("Hartigan-Wong"), trace = FALSE)

data_$cluster <- as.factor(k$cluster)
# data_$cluster <- as.factor(gmm$classification)
# data_$is_elite <- as.factor(data_$is_elite)

write.table(data_, file.path(project_dir, outputs_dir, "nucleotide_transformer/clustering",
  # paste0(timestamp, "summarised_by_enhancer_clusteredGMM_modelVEV_3cl_InclMedian.txt")),
  paste0(timestamp, "summarised_by_enhancer_clusteredkmeans_3cl_InclMedianMaxDist.txt")),
   quote= FALSE, sep = "\t")

# Saving models
saveRDS(k, file.path(project_dir, outputs_dir, "nucleotide_transformer/clustering",
  # paste0(timestamp, "GMM_modelVEV_3cl_InclMedian.rds")))
  paste0(timestamp, "model_kmeans_3cl_InclMedianMaxDist.rds")))

ggplot(data_, aes(x=median_EP_dist, y = enh_count, color = cluster)) +
  geom_point() +
  scale_color_manual(values = brewer.pal(4, "Set1"))

# DENSITY
ggplot(data_, aes(x = median_EP_dist, fill = cluster)) +
  geom_density(alpha = 0.5) 

ggsave(file.path(project_dir, outputs_dir, "nucleotide_transformer/clustering",
  paste0(timestamp, "density_kmeans_3cl_InclMedianMaxDist.png") ))

# PAIRS
ggpairs(select(data_, -gh_pos), aes(color=cluster, alpha = 0.3)) +
  theme_bw() +
  scale_color_manual(values = brewer.pal(8,"Set1"))

ggsave(file.path(project_dir, outputs_dir, "nucleotide_transformer/clustering",
  paste0(timestamp, "pairs_plots_kmeans_3cl_InclMedianMaxDist.png")),
  width = 15, height= 15)


lr_cluster <- ifelse(median(filter(data_, cluster ==1)$median_EP_dist) > median(filter(data_, cluster ==2)$median_EP_dist), 1,2 )
sr_cluster <- ifelse(lr_cluster == 2, 1, 2)

# cutoff_long <- filter(data_, cluster == 1) %>% pull(median_EP_dist) %>% median()
cutoff_long <- filter(data_, cluster == lr_cluster) %>% pull(median_EP_dist) %>% quantile(probs = 0.7)

# cutoff_short <- filter(data_, cluster == 2) %>% pull(median_EP_dist) %>% median()
cutoff_short <- filter(data_, cluster == sr_cluster) %>% pull(median_EP_dist) %>% quantile(probs = 0.15)

ggplot(data_, aes(x = median_EP_dist, fill = cluster)) +
  geom_density(alpha = 0.5) +
  # geom_vline(xintercept = cutoff_long, linetype="dashed", color = "purple", size=0.9) +
  # geom_vline(xintercept = cutoff_short, linetype="dashed", color = "purple", size=0.9) +
  xlim(0,10**6) +
  ylim(0, 0.00005) +
  scale_color_manual(values = brewer.pal(2, "Set1"))

# Defining long and short ranges
# nucleotide_transformer_gr <- data_ %>%
#     mutate(class = case_when(median_EP_dist > cutoff_long & cluster == lr_cluster ~ "long-range",
#         median_EP_dist < cutoff_short & cluster == sr_cluster ~ "short-range",
#         TRUE ~ "not-defined"))

nucleotide_transformer_gr <- data_ %>%
    mutate(class = case_when( cluster == 3 ~ "long-range",
         cluster == 1 ~ "short-range",
        TRUE ~ "not-defined"))

# write.table(nucleotide_transformer_gr, file.path(project_dir, outputs_dir, "nucleotide_transformer/enhancers_annotation_all_gr.txt"), quote=FALSE,row.names=FALSE, sep = "\t")


nucleotide_transformer_gr <- grangesFromDataframe(nucleotide_transformer_gr, colName = gh_pos)

# Filtering out intermediates 
nucleotide_transformer_gr <- nucleotide_transformer_gr[nucleotide_transformer_gr$class != "not-defined"]

# Filtering out width less than 300
min_enh_width <- 300
table(width(nucleotide_transformer_gr) > min_enh_width)
table(width(filter(nucleotide_transformer_gr, class != "not-defined")) > min_enh_width)
table(width(filter(nucleotide_transformer_gr, class == "long-range")) > min_enh_width)
table(width(filter(nucleotide_transformer_gr, class == "short-range")) > min_enh_width)

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

# Filtering width for nucleotide transformer
short_long_gr <- nucleotide_transformer_gr[width(nucleotide_transformer_gr)>=300]

short_long_seq <- as.data.frame(getSeq(BSgenome.Hsapiens.UCSC.hg19, short_long_gr[short_long_gr$class != "not-defined"]))
short_long_seq$label <- ifelse(short_long_gr[short_long_gr$class != "not-defined"]$class == "long-range", "long", "short")
write.table(short_long_seq %>%  select(label, sequence = "x"),
  file.path(project_dir, outputs_dir, "nucleotide_transformer/enh_sequences_short_long_labeled_width300.hg19.csv"),
  quote=FALSE, row.names = FALSE, sep = ",")


