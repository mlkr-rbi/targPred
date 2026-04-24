args = commandArgs(trailingOnly=TRUE)

library(GenomicRanges)
library(dplyr)
library(stringr)
library(stats)
library(parallel)


print("Helloworld")

source("src/empP/empP_utils.R")



# Number of cores
num_cores <- 10

cl <- makeCluster(num_cores)


# Max distance 6.4 Mb
max_dist <- 6400000

# Inputs

## enhancers tpm pooled gr
fantom5.enhancers.tpm.gr.fn <- "data/inputs/fantom5.enhancers.tpm.pooled.gr.rds"
## promoter cage tpm df
promoter_cage_tpm_df.fn <- "data/inputs/promoter_cage_tpm_df_fixed_pooled.rds"
## transcripts
transcripts.fn <- "data/inputs/TSSs.perTranscript.gr.rds"

# Outputs



#Log
log_fn <- file.path("data/processed/logs", "empP_log.txt")

# Loading data
## storage
# dir_precalced = "/storage/PROJECTS/targPred/data/processed"
# fn_out = "/home/kprskalo/targPred/calculating_empP"
# fn_out <- "/home/kprskalo/targPred/calculating_empP/empP_et_all_files_denovo"
# fn_out <- args[2]

chr <- args[1]

## output dir
fn_out <- "data/processed/empP"
fn_out <- file.path(fn_out, chr)

# Create output directory per chromosome if it doesn't exist
if(!dir.exists(fn_out)){
  dir.create(fn_out, recursive = TRUE)
}

#TEST
# fn <- paste0("empP_CAGEenhPermissive_hg19Ens_6.4Mb_",chr,".txt")
# fn <- paste0("02-test_empP_CAGEenhPermissive_hg19Ens_6.4Mb_",chr,".txt")

#TEST
# Switching this to false will load the fn file and continue from the first NA 
create_new_ep = TRUE

# Load data
## These are pooled by cell/tissue type
fantom5.enhancers.tpm.gr <- readRDS(file.path(fantom5.enhancers.tpm.gr.fn))
promoter_cage_tpm_df <- readRDS(file.path(promoter_cage_tpm_df.fn))
transcripts <- readRDS(file.path(transcripts.fn))



print(paste("Processing chromosome", chr))

# Filter out promoters which do not have any expression in the promoter_cage_tpm_df
promoter_cage_tpm_df <- promoter_cage_tpm_df[rowSums(promoter_cage_tpm_df) > 0, ]

empP_df <- CreateEnhPromCombinationsPerChromosome(fantom5.enh.gr = fantom5.enhancers.tpm.gr,
  TSS.gr = transcripts,
  promoter.df = promoter_cage_tpm_df,
  chr = chr,
  max_distance = max_dist, 
  metrics_columns = c("empP_median_oneTailed",
    "empP_median_twoTailed",
    "E_active",
    "eff_size"),
  create_new = create_new_ep,
  filename = file.path(fn_out, fn)
)


empP_df <- CreateEnhPromCombinationsPerChromosome(fantom5.enh.gr = fantom5.enhancers.tpm.gr,
  TSS.gr = transcripts,
  promoter.df = promoter_cage_tpm_df,
  chr = chr,
  max_distance = max_dist, 
  metrics_columns = c("empP_median_oneTailed",
    "empP_median_twoTailed",
    "E_active",
    "eff_size"),
  create_new = create_new_ep,
  filename = file.path(fn_out, fn)
)

# Log the number of enhancer-promoter pairs tested for this chromosome
write(sprintf("Number of enhancer-promoter pairs tested for chromosome %s: %d", chr, nrow(empP_df)), file = log_fn, append = TRUE)

if(create_new_ep){
  empP_df$empP_median_oneTailed <- NA
  empP_df$empP_median_twoTailed <- NA
  # empP_df$empP_ttest_oneTailed <- NA
  # empP_df$empP_ttest_twoTailed <- NA
  # empP_df$empP_wtest_oneTailed <- NA
  # empP_df$empP_wtest_twoTailed <- NA
  empP_df$E_active <- NA
  empP_df$eff_size <- NA
  # empP_df$spearman_r <- NA
  # empP_df$spearman_p <- NA
}

empP_df <- grangesFromDataframe(empP_df, colName = enh_pos)

print(head(empP_df))

i_start <- which(is.na(empP_df$empP_median_oneTailed))[1]
len_ <- length(unique(empP_df$gene))

iteratingEnhPromPairs <- function(id){

  # print(id)
  i <- as.numeric(id)

  # print the progress every 1000 iterations and also the percentage of completion
  if(i %% 1000 == 0){
    print(paste("Processing enhancer-promoter pair", i, "out of", nrow(empP_df)))
    print(paste("Percentage completed:", round(i/nrow(empP_df)*100, 2), "%"))
  }

  enh_expression <- EnhExpression(fantom5.enh.gr = fantom5.enhancers.tpm.gr, enh.gr = empP_df[i])
  
  # print(enh_expression)
  expression_pair <- data.frame(
      enh_expression = enh_expression,
      prom = promoter_cage_tpm_df[rownames(promoter_cage_tpm_df) == empP_df[i]$gene])

    emp_out <- calcEmpP_alt(expression_pair, iter = 10000)

  
  # empP_df$chr <- chr
  empP_df[i]$E_active <- emp_out[1]
  empP_df[i]$empP_median_oneTailed <- emp_out[2]
  empP_df[i]$empP_median_twoTailed <- emp_out[3]
  # empP_df[i]$empP_ttest_oneTailed <- emp_out[4]
  # empP_df[i]$empP_ttest_twoTailed <- emp_out[5]
  # empP_df[i]$empP_wtest_oneTailed <- emp_out[6]
  # empP_df[i]$empP_wtest_twoTailed <- emp_out[7]
  empP_df[i]$eff_size <- EffectSizeOfEnhancerActivityCohens(expression_pair)

  # spearman_cor <- cor.test(x =expression_pair$enh_expression,
  #   y = expression_pair$prom,
  #   method = "spearman")
  
  # empP_df[i]$spearman_r <- spearman_cor$estimate
  # empP_df[i]$spearman_p <- spearman_cor$p.value
  
  out_table <- as.data.frame(mcols(empP_df[i]))
  fn <-  paste0(out_table$enh_pos, "_", out_table$gene,".txt")
  # print(paste("saving the: ",fn))

  
  write.table(out_table,row.names=FALSE, quote=FALSE, file.path(fn_out, fn))
}

mclapply(as.list(1:length(empP_df)), mc.cores = num_cores, FUN = iteratingEnhPromPairs)

stopCluster(cl)