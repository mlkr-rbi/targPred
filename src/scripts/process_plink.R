#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)



print(args)
tb <- read.table(paste0("out/",args,".ld"), header = T)

tb_out <-
  tb %>% 
  group_by(SNP_A) %>%
  summarize(end = max(BP_B), start = min(BP_B))

write.table(tb_out, paste0("out/",args,"_LD.txt"), col.names = FALSE)
