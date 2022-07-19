library(tidyverse)
library(data.table)

# Path from working dir to LDproxy output dir
#LDout_dir <- "data/LDproxy"

gwas_snps_LDblocks <- data.frame()

for(chr in list.dirs()[-1]){
  message(chr)
  
  # Path of each queried SNP
  all_paths <- list.files(path = chr, pattern = "*.txt", full.names = T, recursive = T)
  
  # File names stored as a list
  all_filenames <- all_paths %>%
    basename() %>%
    as.list()
  
  # Read all output files into a single list of data frames
  all_content <- all_paths %>% 
    lapply(read.table,
           header = TRUE,
           sep = "\t",
           encoding = "UTF-8")
  
  # Append filenames to corresponding df in list of dataframes (filenames = queried SNP for each df)
  all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)
  
  # Make a single dataframe from list of dataframes
  all_result <- rbindlist(all_lists, fill = T)
  
  # Separate filename to remove .txt extension, remove unnecessary columns, separate coord column into chromosome number and loc (locus coordinates on chromosome), keep only SNPs correlated with the query snp with R2 > 0.8 
  all_results_0.6 <- separate(all_result, V1, c('query_snp', 'rest'), "[.]") %>%
    select(-rest, -Function) %>%
    separate(Coord, c('chromosome', 'loc'), ':') %>%
    transform(loc = as.numeric(loc)) %>%
    as.data.frame(all_results) %>% 
    filter(R2 >= 0.6)
  
  
  all_results_0.6_chr <- select(all_results_0.6, query_snp, chromosome) %>%
    unique()
  
  
  # Get LD Block range
  # Get table with start and end coordinates for LD block of queried SNP
  range_table <- all_results_0.6 %>% 
    group_by(query_snp) %>% 
    summarize(max_loc = max(loc), min_loc = min(loc)) %>%
    left_join(all_results_0.6_chr)
  
  gwas_snps_LDblocks <- rbind(gwas_snps_LDblocks, range_table)
  
}

saveRDS(gwas_snps_LDblocks, file="gwas_snps_LDblocks.rds")
