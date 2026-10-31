library(org.Hs.eg.db)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(tidyverse)
# These are needed for EPplots.R
library(stringr)
library(dplyr)
library(GenomicRanges)
library(RColorBrewer)
#library(beanplot)
library(biomaRt)
library(parallel)

# Manually set the path to utils.R
source(file.path("src/core", "utils.R"))
# Cached gene/ideogram tracks and Ensembl retry helpers
source(file.path("src/core", "gviz_cache.R"))

directory_type <- ""
paths <- set_paths("src/core/paths.txt")
core_dir <- paths$core_dir
work_dir <- paths$work_dir
data_dir <- paths$data_dir
gwas_year <- paths$gwas_year
byGRBlist_name <- paths$byGRBlist_name

# Paths derived from already loaded objects
cne_density_file <- file.path(
  data_dir, "inputs", "HCNE_density_hg19_canFam3_98pc_50col.wig.gz"
)
directionality_file <- file.path(
  data_dir, "inputs", "tad_tracks", "directionality_homer_neural.bigWig"
)
gviz_output_suffix <- "_Gviz.png"

# use set_paths using paths.txt
# set_paths(path_file)

# Object that contains all GRBs, genes, enhancers and empP values for enh-transcript
byGRB_list <- readRDS(file.path(data_dir,"processed",byGRBlist_name))
loci_in_GRBs.df <- read.csv(file.path(data_dir,"processed","loci_in_GRBs.df.csv"), header = TRUE, stringsAsFactors = FALSE)
# loci_in_GRBs.gr <- readRDS(file.path(data_dir,"processed","loci_in_GRBs.gr.rds"))
snps_in_GRBs.gr <- readRDS(file.path(data_dir,"processed","snps_in_GRBs.gr.rds"))
# Target genes for CNE density GRBs 
target_genes <- 
  read.csv(file.path(data_dir,"inputs","hg19.canFam3.tsv"), sep = "\t", header = TRUE)  
# Building target gene genomic ranges
# enhancers and promoters are pooled data
TSSs.perTranscript.gr <- readRDS(file.path(data_dir, "inputs","TSSs.perTranscript.gr.rds"))
fantom5.enhancers.tpm.gr <- readRDS(file.path(data_dir, "inputs","fantom5.enhancers.tpm.pooled.gr.rds"))
promoter_cage_tpm_df <- readRDS(file.path(data_dir, "inputs","promoter_cage_tpm_df_fixed_pooled.rds"))
# GRBS
grbs <- readRDS(file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))
# list with gene ids (Anja)
genes.hg19.ensembl <- readRDS(file.path(data_dir, "inputs","genes.hg19.ensembl.rds"))

# gwas_rds_path <- file.path(data_dir, "inputs", paste0("gwas_catalog_",gwas_year), "gwas_associations.rds")
# message("[diagnostic] R working dir: ", getwd())
# message("[diagnostic] gwas_rds_path: ", gwas_rds_path)
# message("[diagnostic] file.exists: ", file.exists(gwas_rds_path))
# message("[diagnostic] list parent dir: ")
# print(list.files(dirname(gwas_rds_path), all.files = TRUE))
gwas_associations <- read.csv(file.path(data_dir, "outputs", "gwas_associations.csv"), header = TRUE, stringsAsFactors = FALSE)

col_vec <- c("darkgoldenrod1", "darkorange1", "red")

gviz_dir <- file.path( "data/plots/gviz_plots")

# Check if gviz_dir exists, if not create it
if(!file.exists(gviz_dir)){
  dir.create(gviz_dir)
}

# Filter NA for MAPPED_TRAIT in loci_in_GRBs.df
loci_in_GRBs.df <- loci_in_GRBs.df[!is.na(loci_in_GRBs.df$MAPPED_TRAIT),]

# Create loci_in_GRB.gr from df
loci_in_GRBs.gr <- grangesFromDataframe(loci_in_GRBs.df, colName = LDblock, keep_metadata = TRUE)

# Adding a column to grbs that correspond to targets
# actual plotting, saved in gviz_dir
loci_in_GRBs.gr$name <- loci_in_GRBs.gr$SNPS # need to have a $name column for gviz plotting



# Gene tracks come from the GRCh37 Ensembl archive (do NOT use the default, hg38, mart).
# The Mart connection, gene tracks and ideograms are cached under data/cache/gviz;
# the connection is only opened on a cache miss (see gviz_cache.R).
options(targpred.gviz_cache_dir = file.path(data_dir, "cache", "gviz"))
extend_biomart_timeout(180)
mart <- NULL


# Looping through traits and plotting Gviz for each trait separately, including only loci with the trait and SNPs with the trait
traits <- unique(loci_in_GRBs.df$MAPPED_TRAIT) %>% sort()

num_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = 1))

# Each trait returns one status row per GRB panel: written, skipped (the helper
# found no usable loci/enhancers) or error.
status_list <- mclapply(traits, function(trait) {
    use_local_biomart_cache()
    trait_dir <- file.path(gviz_dir,str_replace_all(trait," |/","_"))
    # Creates a dir if it does not exist
    if(!file.exists(trait_dir)){dir.create(trait_dir)}

    # to be done on each locus
    message(trait)
    # Including only loci in trait
    loci_in_GRB_trait <- loci_in_GRBs.gr[loci_in_GRBs.gr$MAPPED_TRAIT == trait,]
    # Remoivng duplicated SNPS per trait/locus
    loci_in_GRB_trait <- loci_in_GRB_trait[!duplicated(loci_in_GRB_trait$SNPS),]

    grbs_to_loop <-
        which( names(byGRB_list) %in% str_replace_all(loci_in_GRB_trait$grb.coords, pattern = ":|-", replacement = "_"))

    status <- data.frame(trait=character(), grb=character(), status=character(), detail=character())
    #for (i in (1:length(byGRB_list))) {
    for (i in grbs_to_loop) {
        print(paste("Plotting Gviz for", names(byGRB_list)[i], sep=" "))
        plot_dir <- paste0(trait_dir,"/",str_replace_all(trait," |/","_"),"_")
        output <- paste0(plot_dir, str_replace_all(names(byGRB_list)[i], ":", "_"), gviz_output_suffix)
        previous_mtime <- file.mtime(output)
        detail <- ""
        outcome <- tryCatch({
            with_retry(function() tryCatch(
                locusGviz(i,
                          dir=plot_dir,
                          labels=T,
                          byGRB_list=byGRB_list,
                          # Selecting only loci with the trait
                          loci.gr= loci_in_GRB_trait,
                          # Selecting only SNPs with the trait
                          snps.gr=snps_in_GRBs.gr[snps_in_GRBs.gr$SNPS %in% loci_in_GRB_trait$SNPS,],
                          grbs=grbs,
                          mart=mart,
                          trait = trait,
                          cne_density_file = cne_density_file,
                          directionality_file = directionality_file,
                          output_suffix = gviz_output_suffix),
                error=function(e) { while (length(dev.list())) dev.off(); stop(e) }))
            # A PNG left untouched from an earlier run does not count as written
            rewritten <- file.exists(output) && (is.na(previous_mtime) || file.mtime(output) != previous_mtime)
            if (rewritten && file.size(output) > 1000) "written" else "skipped"
        }, error=function(e) {
            detail <<- conditionMessage(e)
            cat("ERROR :", detail, "\n", trait, i, "\n")
            "error"
        })
        status <- rbind(status, data.frame(trait=trait, grb=names(byGRB_list)[i],
                                           status=outcome, detail=detail))
    }
    status
}, mc.cores = num_cores)

# A worker that died entirely returns a try-error (or NULL if it was killed)
# instead of a data frame
worker_failed <- vapply(status_list, function(x) is.null(x) || inherits(x, "try-error"), logical(1))
failed_traits <- traits[worker_failed]
status_df <- bind_rows(status_list[!worker_failed])
n_failed <- length(failed_traits)
status_df <- bind_rows(status_df, data.frame(trait=failed_traits, grb=rep(NA_character_, n_failed),
                                             status=rep("error", n_failed),
                                             detail=rep("worker failed", n_failed)))
status_file <- file.path(gviz_dir, "gviz_status.tsv")
write.table(status_df, status_file, sep="\t", quote=TRUE, row.names=FALSE)
print(table(status_df$status))
if (any(status_df$status == "error")) {
    stop(sum(status_df$status == "error"), " Gviz panels failed; see ", status_file)
}
