library(GenomicRanges)
library(tidyverse)

project_dir <- "/home/katarina/MyProjects/targPred"
inputs_dir <- "data/inputs"
processed_dir <- "data/processed"
outputs_dir <- "/data/outputs"

source(file.path(project_dir, "src/scripts/utils.R"))

# GRBS
grb <- readRDS(file.path(project_dir, processed_dir, "GRBs_with_targets.gr.rds"))


fn <-file.path(project_dir, processed_dir, "genome_browser/modified_GRBs_with_targets.bed")
makeBed(grb, name = grb$targets, filename = fn, includetrackname = TRUE,
        track_name = "Grb modified",
        track_description = "GRBs CNE density modified with kurtosis method")

# Loci in GRB
loci_in_GRBs.gr <- readRDS(file.path(project_dir, processed_dir, "loci_in_GRBs.gr.rds"))
fn <-file.path(project_dir, processed_dir, "genome_browser/loci_in_GRBs.bed")
makeBed(loci_in_GRBs.gr, name = loci_in_GRBs.gr$SNPS, filename = fn, includetrackname = TRUE,
        track_name = "Loci in GRB",
        track_description = "LD blocks in European ancestry for disease associated SNPs")
