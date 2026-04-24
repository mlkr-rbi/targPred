# Crude numbers analysis for targPred paper


# Libraries
library(stringr)
library(dplyr)
library(readr)
library(ggplot2)
library(tibble)
library(ggpubr)
library(viridis)
library(GenomicRanges)

project_dir <- getwd()
inputs_dir <- "data/inputs"
processed_dir <- "data/processed"
outputs_dir <- "data/outputs"


source(file.path(project_dir, "src/core/utils.R"))

paths <- parse_paths_txt("src/core/paths.txt")
empP_dir <- paths$empP_dir
empP_basename <- paths$empP_basename
gwas_year <- paths$gwas_year
LDblock_path <- paths$LDblock_path
byGRBlist_name <- paths$byGRBlist_name

# Define input and output directories
inputs <- "data/analysis/crude_numbers/inputs"
outputs <- "data/analysis/crude_numbers/outputs"

# INPUT FILES:
gwas_fn <- file.path(project_dir, "data/outputs/gwas_associations.csv")
web_table_fn <- file.path(project_dir, "data/outputs/web_table.csv")
TSSs.perTranscript.gr <- readRDS(file.path(inputs_dir, "TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr_fn <- file.path(inputs_dir, "fantom5.enhancers.tpm.pooled.gr.rds")
promoter_cage_tpm_df <- readRDS(file.path(inputs_dir, "promoter_cage_tpm_df_fixed_pooled.rds"))
loci_EP_fn <- file.path(outputs, "loci_EP.csv")


# OUTPUT FILES:
#LOG file
log_fn <- file.path(project_dir, processed_dir, "logs", "crude_numbers.txt")

# - /home/katarina/MyProjects/targPred/data/processed/targPred_exploration_cumulative_enhancer_df.rds
# - /home/katarina/MyProjects/targPred/data/outputs/byGRB_list_twosided.rds
# - /home/katarina/MyProjects/targPred/data/inputs/hg19.canFam3.tsv
# - /home/katarina/MyProjects/targPred/data/processed/empP_files_twosided/
#
# OUTPUT FILES:
# - /home/katarina/MyProjects/targPred/data/outputs/targPred_exploration/nested_pie_plot_white_v1.pdf
# - /home/katarina/MyProjects/targPred/data/outputs/targPred_exploration/number_of_enhancers_cumulative_binned_distance_target_gene_v3.pdf
# - /home/katarina/MyProjects/targPred/data/outputs/targPred_exploration/number_of_enhancers_cumulative_binned_distance_target_gene_test.csv

# LOAD DATA
gwas <- read.csv(gwas_fn)
web_table <- read.csv(web_table_fn)
enhancers <- readRDS(fantom5.enhancers.tpm.gr_fn)

enhancer_transcript_df <- read.csv(loci_EP_fn)

# Colors
GRB_CAT_COLORS <- c("GRB-byst" ="#143159", 
                    "GRB-trgt" = "#4E7DD9",
                    "GWAS-byst" = "#D9A60D",
                    "GWAS-out" = "#275912",
                    "GWAS-trgt" = "#DA3426")

# Basic results calculation
basic_results <- c()
basic_results["GWAS_processed"] <- dim(web_table)[1]
basic_results["Loci_in_GRB"] <- table(web_table$inGRB)["TRUE"]
basic_results["Loci_out_GRB"] <- table(web_table$inGRB)["FALSE"]

# GWAS and GRB annotation match
genes2 <- strsplit(web_table$GWAS_genes, ", ")
genes1 <- strsplit(web_table$GRB_targets, ", ")
result <- sapply(seq_along(genes1), function(i) {
  any(genes1[[i]] %in% genes2[[i]])
})
data_ <- web_table
data_$matches_target_gene <- result
data_ <- data_ %>%
  filter(inGRB == TRUE) %>%
  filter(GRB_targets != "unknown") 
basic_results["GWAS_and_GRB_annotation_match"] <- table(data_$matches_target_gene)["TRUE"]
basic_results["GRB_novel_annotation"] <- table(data_$matches_target_gene)["FALSE"]

# Loci overlapping enhancers

loci_inGRB.gr <- grangesFromDataframe(filter(web_table, inGRB == TRUE), colName = Locus_position)
loci_outGRB.gr <- grangesFromDataframe(filter(web_table, inGRB == "FALSE"), colName = Locus_position)
loci_overlap_enh <- subsetByOverlaps(loci_inGRB.gr, enhancers)$SNP
nonGRBloci_overlap_enh <- subsetByOverlaps(loci_outGRB.gr, enhancers)$SNP

data_ <-  web_table %>% filter(is.element(SNP, loci_overlap_enh))
basic_results["Loci_in_GRB_overlapping_enhancers"] <- dim(data_)[1]
data_ <-  web_table %>% filter(is.element(SNP, nonGRBloci_overlap_enh))
basic_results["Loci_not_in_GRB_overlapping_enhancers"] <- dim(data_)[1]

# Bar plot
df <- data.frame(name = names(basic_results), value = basic_results)
df$name <- factor(df$name, levels = names(basic_results))
ggplot(df, aes(x = name, y = value, fill = name)) +
  geom_bar(stat = "identity") +
  scale_fill_viridis_d() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Count", title = "Basic Results")

# Stacked bar plot
df <- data.frame(cat = factor(c("Number of GWAS loci", "Number of GWAS loci", "Loci overlapping enhancers", "Loci overlapping enhancers",  "Gene annotation", "Gene annotation"), levels = c("Number of GWAS loci", "Loci overlapping enhancers", "Gene annotation")), subcat = c("in GRB", "out GRB", "in GRB", "out GRB", "GWAS and GRB match", "GRB novel"), value = basic_results[c(2,3,6,7,4,5)])
col_vec <- c("out GRB" = "#143159", "in GRB" = "#4E7DD9", "GRB novel" = "#D9A60D", "GWAS and GRB match" = "#DA3426")
ggplot(df, aes(x = cat, fill = subcat, y = value)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  scale_fill_manual(values = col_vec) +
  theme_bw() +
  labs(x = NULL, y = NULL) +
  theme(legend.title = element_blank())

# Nested Pie plot
library(tidyr)
middle_ring <- df %>% filter(cat == "Number of GWAS loci") %>% select(cat, subcat, value)
## inner ring
inner_ring <- df %>% filter(cat == "Loci overlapping enhancers") %>% select(cat, subcat, value)
loci_in_GRB_no_ov <- basic_results["Loci_in_GRB"] - inner_ring$value[inner_ring$subcat == "in GRB"]
loci_out_GRB_no_ov <- basic_results["Loci_out_GRB"] - inner_ring$value[inner_ring$subcat == "out GRB"]
inner_ring <- inner_ring %>%
  rbind(c("Loci overlapping enhancers", "in GRB - no ov", loci_in_GRB_no_ov)) %>%
  rbind(c("Loci overlapping enhancers", "out GRB - no ov", loci_out_GRB_no_ov))
inner_ring$value <- as.numeric(inner_ring$value)
inner_ring$subcat <- factor(inner_ring$subcat, levels = rev(c("out GRB", "out GRB - no ov", "in GRB - no ov", "in GRB")))
## outer ring
outer_ring <- df %>% filter(cat == "Gene annotation") %>% select(cat, subcat, value) 

annotation_exists <- basic_results["GWAS_and_GRB_annotation_match"] + basic_results["GRB_novel_annotation"]
outer_ring <- outer_ring %>%
  rbind(c("Gene annotation", "other", basic_results["Loci_out_GRB"])) %>%
  rbind(c("Gene annotation", "unknown", basic_results["Loci_in_GRB"] - annotation_exists))
outer_ring$value <- as.numeric(outer_ring$value)
outer_ring$subcat <- factor(outer_ring$subcat, levels = c("GWAS and GRB match", "GRB novel", "unknown", "other"))
other_col <- "white"
col_vec <- c("out GRB" = "#143159", "in GRB" = "#4E7DD9", "GRB novel" = "#D9A60D", "GWAS and GRB match" = "#DA3426", "in GRB - no ov" = other_col, "out GRB - no ov" = other_col, "other" = other_col, "unknown" = other_col)
plt <- ggplot() + 
  geom_col(aes(x = 1, y = value, fill = subcat), data = inner_ring ) +
  geom_col(aes(x = 2, y = value, fill = subcat), data = middle_ring) +
  geom_col(aes(x = 3, y = value, fill = subcat), data = outer_ring) +
  xlim(0, 3.5) + 
  labs(x = NULL, y = NULL) + 
  scale_fill_manual(values = col_vec) +
  theme_minimal() +
  theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank(), legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank())
plt <- plt + coord_polar(theta = "y")
# ggsave(height = 8, width = 7, filename = file.path(project_dir, outputs_dir, "targPred_exploration", "nested_pie_plot_white_v1.pdf"))
svg(filename = file.path(outputs, "nested_pie_plot_white_v1.svg"), height = 8, width = 7)
print(plt)
dev.off()

# Ratios to Log
cat(paste0("Ratio of loci that fall in GRBs ", basic_results["Loci_in_GRB"]/(basic_results["Loci_in_GRB"]+basic_results["Loci_out_GRB"])), "\n", file = log_fn, append = TRUE)
# ratio of loci that fall in GRBs and overlap enhancers
cat(paste0("Ratio of loci that fall in GRBs and overlap enhancers ", basic_results["Loci_in_GRB_overlapping_enhancers"]/basic_results["Loci_in_GRB"]), "\n", file = log_fn, append = TRUE)
cat(paste0("Ratio of novel gene annotation ", basic_results["GRB_novel_annotation"]/basic_results["Loci_in_GRB"]), "\n", file = log_fn, append = TRUE)

# get transcripts with expression across promoters == 0
null_promoters <- rownames(promoter_cage_tpm_df[rowSums(promoter_cage_tpm_df) == 0, ])

# split gene:transcript names to keep only transcript names
null_promoters <- sapply(strsplit(null_promoters, split = ":"), function(x) x[2])

# extract list of genes with multiple transcripts
multi_transcript_genes <- enhancer_transcript_df %>%
  select(ensembl, transcript) %>%
  unique() %>%
  group_by(ensembl) %>%
  count() %>%
  filter(n > 1) %>%
  pull(ensembl)

# extact transcripts from multi-transcript genes to keep only transcripts with the largest number of cm_enh_num_4gene
multi_transcripts_to_filter <- enhancer_transcript_df %>%
  filter(ensembl %in% multi_transcript_genes) %>%
  group_by(ensembl) %>%
  filter(cm_enh_num_4gene == max(cm_enh_num_4gene)) %>%
  ungroup() %>%
  select(transcript) %>%
  pull()


# Number of enhancers per gene
give.n <- function(x){ return(c(y = -0.5, label = length(x))) }
enhtrans_sub <- enhancer_transcript_df %>%
  # Using FDR pajd
  # filter(p_adj < 0.05) %>%
  select(ensembl, symbol, transcript, ranges_binned, target_gene, cm_enh_num_4gene) %>%
  unique() %>%
  filter(!is.na(ranges_binned)) %>%  
  # filter out cm_enh_num_4gene == 0 (if a transcript has no significant EP pairs, the cumulative count should be 0)
  filter(cm_enh_num_4gene > 0) %>%
  # filter out multi-transcripts to keep only the one with the largest number of cm_enh_num_4gene
  # filter(!(transcript %in% multi_transcripts_to_filter))
  # filter out null_promoters
  filter(!(transcript %in% null_promoters))


p <- ggplot(enhtrans_sub, aes(x = target_gene, y = cm_enh_num_4gene, fill = target_gene)) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  facet_wrap(~ranges_binned, scales = "free", nrow = 1, ncol = 5) +
  coord_cartesian(ylim = c(-0.8, 30)) +
  scale_fill_manual(values = c("#D9A60D","#4E7DD9")) +
  stat_summary(fun.data = give.n, geom = "text") +
  labs(x = "", fill = "Target gene", y = "Number of enhancers") +
  theme(plot.background = element_rect(fill = "white"), axis.ticks = element_blank(), strip.text.y = element_text(size = 12), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

# Open a PDF device to save the plot
svg(file.path(outputs, "number_of_enhancers_cumulative_binned_distance_target_gene_v7.svg"), height = 6, width = 9)
print(p)
dev.off()

# Wilcoxon test
wilcox_summary_tbl <- data.frame(ranges_binned = character(), median_target_gene = numeric(), median_bystander_gene = numeric(), W = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
for(i in enhancer_transcript_df$ranges_binned %>% unique() %>% sort){
  x = filter(enhtrans_sub, ranges_binned == i & target_gene == "No")$cm_enh_num_4gene
  y = filter(enhtrans_sub, ranges_binned == i & target_gene == "Yes")$cm_enh_num_4gene
  wt <- wilcox.test(x, y, alternative = "less")
  wilcox_summary_tbl <- rbind(wilcox_summary_tbl, data.frame(ranges_binned = i, median_No = median(x, na.rm = TRUE), median_Yes = median(y, na.rm = TRUE), W = unname(wt$statistic), p_value = wt$p.value, stringsAsFactors = FALSE))
}
write_csv(wilcox_summary_tbl, file =  file.path(outputs, "number_of_enhancers_cumulative_binned_distance_target_gene_test.csv"))

# t test
for(i in enhancer_transcript_df$ranges_binned %>% unique() %>% sort){
  x = filter(enhtrans_sub, ranges_binned == i & target_gene == "No")$cm_enh_num_4gene
  y = filter(enhtrans_sub, ranges_binned == i & target_gene == "Yes")$cm_enh_num_4gene
  print(t.test(x,y, alternative = "less"))
}

# Session Info
sessionInfo()
