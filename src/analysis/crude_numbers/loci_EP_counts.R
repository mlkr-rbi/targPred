# targPred EP exploration - R script version
# INPUT FILES:
# - byGRB_list_twosided.rds
# - enh_trans_ranges.rds
# - hg19.canFam3.tsv
# OUTPUT FILES:
# - by_GRB_list_twosided_df.rds
# - targPred_exploration_cumulative_enhancer_df.rds
# - plots/enh_num_4gene_target_vs_nontarget.png
# - targPred_exploration/number_of_enhancers_cumulative_binned_distance_target_gene_v2.png

# Libraries
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(GenomicRanges)

project_dir <- getwd()
inputs_dir <- "data/inputs"
processed_dir <- "data/processed"
outputs_dir <- "data/outputs"


source(file.path(project_dir, "src/core/utils.R"))

paths <- parse_paths_txt("src/core/paths.txt")
empP_dir <- paths$empP_dir
empP_basename <- paths$empP_basename
byGRBlist_name <- paths$byGRBlist_name

# Define input and output directories
inputs <- "data/analysis/crude_numbers/inputs"
outputs <- "data/analysis/crude_numbers/outputs"

# INPUT FILES:
## From main pipeline:
byGRB_list_fn <- file.path(project_dir, processed_dir, byGRBlist_name)
grb_with_targets_fn <- file.path(project_dir, processed_dir, "GRBs_with_targets.gr.rds")


# OUTPUT FILES:
#LOG file
log_fn <- file.path(project_dir, processed_dir, "logs", "crude_numbers.txt")

# Enhancer-transcript distance file
EP_dist_fn <- file.path(project_dir, outputs, "EP_dist_df.csv")
# Main table with enhancer-transcript pairs and their empP values and number of interactions
loci_EP_fn <- file.path(project_dir, outputs, "loci_EP.csv")


# Load data
byGRB_list <- readRDS(byGRB_list_fn)
# trg_genes <- read.table(file.path(project_dir, inputs_dir , "hg19.canFam3.tsv"), header = TRUE)
grb <- readRDS(grb_with_targets_fn)


# Globals
P_VAL <- 0.05

GRB_CAT_COLORS <- c("GRB-byst" ="#143159", 
                    "GRB-trgt" = "#4E7DD9",
                    "GWAS-byst" = "#D9A60D",
                    "GWAS-out" = "#275912",
                    "GWAS-trgt" = "#DA3426")

# Helper functions
namedVecToDf <- function(named_vector, column_name = "col1"){
  df <- as.data.frame(named_vector) %>%
    dplyr::select({{column_name}} := Var1, Freq)
  return(df)
}

calcHistogramDataMaxBin <- function(data_, max_bin){
  data_tb <- table(data_) %>% as.data.frame()
  total_ <- data_tb[(max_bin+1):dim(data_tb)[1], "Freq"] %>% sum()
  res_ <- data_tb[0:max_bin,] %>%
    mutate(data_= as.character(data_))
  res_ <- rbind(res_,c("rest", total_)) %>%
    select(x=data_, Freq)
  res_$x <- factor(res_$x, levels = res_$x)
  return(res_)
}

give.n <- function(x){
  return(c(y = -0.5, label = length(x))) 
}

# Load empP files to get transcript locations and append to the enhancer_transcript_df table
EP_dist_df <- data.frame()

for(i in c(1:22, "X")){
  empP_file <- file.path(project_dir, processed_dir, empP_dir, paste0(empP_basename, "chr", i, ".txt"))
  empP_chr <- read.table(empP_file, header = TRUE, sep = " ")

  # Assuming your data frame is called df
  empP_chr <- empP_chr %>%
    select( enh_pos, gene, transcript_position) %>%
    unique() %>%
    # mutate(transcript_position = floor((start + end) / 2)) %>%
    # select(transcript, chr, transcript_position)
  # Parse enhancer start and end from enh_pos
  mutate(
    enh_start = as.numeric(str_extract(enh_pos, "(?<=:)[0-9]+")),
    enh_end   = as.numeric(str_extract(enh_pos, "(?<=-)[0-9]+")),
    enh_mid   = floor((enh_start + enh_end) / 2),
    EP_dist   = abs(enh_mid - as.numeric(transcript_position))
  )
  # select(chr, enh_pos, gene, transcript_position, EP_dist, everything()) %>%
  # distinct(enh_pos, gene, transcript_position, .keep_all = TRUE)

  EP_dist_df <- rbind(EP_dist_df, select(empP_chr, enh_pos, transcript_position, gene, EP_dist))
  
  print(paste("Finished processing chromosome", i))
  # head and dim of enhancer_transcript_df after processing

}


#LOG number of rows in EP_dist_df
cat("Number of rows in EP_dist_df which represents the total number of EP pairs tested in empP files:", nrow(EP_dist_df), "\n", file = log_fn, append = TRUE)

# write EP_dist_df to a csv file for future use
write.csv(EP_dist_df, file = EP_dist_fn, row.names = FALSE)


enhancer_transcript_df <- data.frame()
empty_grbs <- c()

for(i in 1:length(byGRB_list)){
  message(paste(i,". Binding the grb: ",names(byGRB_list)[i], sep = ""))
  working_el <- byGRB_list[[i]]
  if(setequal(dim(byGRB_list[[i]][[3]]) ,c(0,0))){
    message(paste(i,". does not have enhancers", sep=""))
    empty_grbs <- append(empty_grbs, names(byGRB_list)[i])
    next
  }
  working_df <- 
    data.frame(
      grb = names(byGRB_list)[i],
      gene_trans_enh = c(outer(colnames(working_el[[3]]),rownames(working_el[[3]]), FUN=paste, sep="_")),
      emp = as.vector(t(as.matrix(working_el[[3]])))
    ) %>%
    dplyr::mutate(gene_trans=str_split(gene_trans_enh, pattern="_", simplify=TRUE)[,1],
                  enh=str_split(gene_trans_enh, pattern="_", simplify=TRUE)[,2]) %>%
    dplyr:: mutate(ensembl=str_split(gene_trans, pattern=":", simplify=TRUE)[,1]) %>%
    dplyr::mutate(transcript=str_split(gene_trans, pattern=":", simplify=TRUE)[,2]) %>%
    left_join(working_el[[1]], by=c("ensembl"="ens")) %>%
    dplyr::select(-gene_trans_enh) %>%
    # Filter out GWAS_out genes using category column
    filter(category != "GWAS-out")

  working_df <- working_df %>%
    group_by(grb) %>%
    mutate(p_adj = p.adjust(emp, method = "fdr")) %>%
    ungroup()

  if(any(working_df$p_adj < P_VAL)){
    working_df <- working_df %>%
      left_join(as.data.frame(base::table(filter(working_df, p_adj < P_VAL)$enh)), 
                by = c("enh" = "Var1")) %>%
      dplyr::rename(gene_num_4enh=Freq) %>%
      left_join(as.data.frame(base::table(filter(working_df, p_adj < P_VAL)$transcript)), 
                by = c("transcript" = "Var1")) %>%
      dplyr::rename(enh_num_4gene=Freq)
  } else {
    working_df$gene_num_4enh <- NA
    working_df$enh_num_4gene <- NA
  }
  enhancer_transcript_df <- rbind(enhancer_transcript_df, working_df)
}


# enhancer_transcript_df <-
#   enhancer_transcript_df %>%
#   dplyr::left_join(EP_dist_df, by = c("enh"= "enh_pos","transcript"))
# enh_midpoint <- floor((as.numeric(str_split(enhancer_transcript_df$enh, pattern = ":|-", simplify = TRUE)[,2]) + as.numeric(str_split(enhancer_transcript_df$enh, pattern = ":|-", simplify = TRUE)[,3]))/2)

enhancer_transcript_df <- enhancer_transcript_df %>%
  mutate(gene_num_4enh = ifelse(is.na(gene_num_4enh), 0, gene_num_4enh),
          enh_num_4gene = ifelse(is.na(enh_num_4gene), 0, enh_num_4gene))



enhancer_transcript_df <- enhancer_transcript_df %>%
left_join(select(EP_dist_df, enh_pos, gene, EP_dist), by = c("enh" = "enh_pos", "gene_trans" = "gene")) 

#LOG dim of enhancer_transcript_df after joining with EP_dist_df
cat("Dimensions of enhancer_transcript_df after joining with EP_dist_df:", dim(enhancer_transcript_df), "\n", file = log_fn, append = TRUE)


enhancer_transcript_df <-
  enhancer_transcript_df %>%
  dplyr::mutate(ranges_binned = case_when(
    EP_dist <= 1000 ~ "1-10^3",
    EP_dist > 1000 & EP_dist <= 10000 ~ "10^3-10^4",
    EP_dist > 10000 & EP_dist <= 100000 ~ "10^4-10^5",
    EP_dist > 100000 & EP_dist <= 1000000 ~ "10^5-10^6",
    EP_dist > 1000000 ~ "10^6-"))

# EP pair count binned
EP_pair_count <- 
  enhancer_transcript_df %>%
  filter(p_adj < P_VAL) %>%
  group_by(transcript, ranges_binned) %>%
  count(transcript) %>%
  select(enh_num_4gene_binned = n, everything())

enhancer_transcript_df <- left_join	(enhancer_transcript_df, EP_pair_count)


# Adding GRB target gene annotation 
grb_targets <- str_split(grb$ENSEMBLGeneID, pattern = ",", simplify = TRUE)
# Exclude empty target gene lists
grb_targets <- grb_targets[grb_targets != " "]
# Exlcude NA values
grb_targets <- grb_targets[!is.na(grb_targets)]

enhancer_transcript_df <- enhancer_transcript_df %>%
  dplyr::mutate(target_gene =
                  case_when(is.element(ensembl, grb_targets) ~ "Yes",
                  !is.element(ensembl, grb_targets) ~ "No"))

# Cumulative enhancer counts for each transcript and binned distance
EP_pair_count <- 
  enhancer_transcript_df %>%
  filter(p_adj < P_VAL) %>%
  select(transcript, ranges_binned, enh_num_4gene_binned) %>%
  unique() %>% 
  arrange(ranges_binned)%>%
  group_by(transcript) %>%
  mutate("cm_enh_num_4gene" = cumsum(enh_num_4gene_binned))

enhancer_transcript_df <- 
  left_join(enhancer_transcript_df, EP_pair_count)

# Replace NA values in cm_enh_num_4gene with 0 (if a transcript has no significant EP pairs, the cumulative count should be 0)
enhancer_transcript_df <- enhancer_transcript_df %>%
  mutate(cm_enh_num_4gene = ifelse(is.na(cm_enh_num_4gene), 0, cm_enh_num_4gene))

# An example of cumulative enhancer counts for a single transcript across binned distance ranges
example <- enhancer_transcript_df %>%
  filter(p_adj < 0.05 & transcript =="ENST00000308873") %>%
  select(symbol, transcript, ranges_binned, enh_num_4gene_binned, cm_enh_num_4gene) %>%
  unique() %>%
  arrange(ranges_binned)

#LOG example of cumulative enhancer counts for a single transcript
cat("Example of cumulative enhancer counts for transcript ENST00000308873:\n", file = log_fn, append = TRUE)
write.table(example, file = log_fn, append = TRUE, row.names = FALSE, col.names = TRUE, sep = "\t")

# Save the final
# saveRDS(enhancer_transcript_df, file = loci_EP)
write.csv(enhancer_transcript_df, file = loci_EP_fn, row.names = FALSE)
