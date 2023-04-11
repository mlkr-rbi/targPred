## Script for processing genome shuffling outputs
library(GenomicRanges)
library(tidyverse)
library(plyranges)

chr_vec <- paste0("chr", c(1:22)) %>%
  append("chrX")
# file_ <- "/home/katarina/MyProjects/targPred/data/random/shuffled_beds/shuffle_491.bed"
file_ <- "/home/katarina/MyProjects/GRB/GRB_thresholds_beds/hg19.canFam3.98.50.processed.threshold.grbs.bed"

ep_dist_vec <- c()

shuffled_grb <- read.table(file = file_)
shuffled_grb <- shuffled_grb[-1,]

shuffled_grb.gr <- GRanges(seqnames = shuffled_grb$V1,
                           ranges = IRanges(start = as.numeric(shuffled_grb$V2) +1 ,
                                            end = as.numeric(shuffled_grb$V3) +1))

# transcripts <- readRDS(file.path(project_dir, processed_dir,  "TSSs.perTranscript.gr.rds"))


# shuffled_grb.gr <- readRDS("/home/katarina/MyProjects/GRB/GRB_thresholds_beds/")

ep_ranges_all <- GRanges()

# for(chr in "chr1"){
for(chr in chr_vec){
  
  # Reading all ep contacts per chromosome
  empP_gr <- readRDS(file = file.path(project_dir, outputs_dir, "empP_EP_distance", paste0(chr, "empP_CAGEenhPermissive_hg19Ens_6.4Mb_with_distances")))

  # Filtering significant ones
  empP_gr <- empP_gr[empP_gr$empP < 0.001]
  
  
  # Filtering enhancers that overlap shuffled grbs
  empP_gr <- subsetByOverlaps(empP_gr, shuffled_grb.gr, type = "within")
  
  # Filtering transcripts that overlap shuffled grbs
  # transcript_ov <- subsetByOverlaps(transcripts, shuffled_grb.gr)
  print(format(Sys.time(), "%X"))
  # empP_gr <- empP_gr[str_detect(empP_gr$gene, paste(transcript_ov$ensembl_transcript_id, collapse = '|'))]
  # empP_gr <- plyranges::filter(empP_gr, str_detect(gene, paste(transcript_ov$ensembl_transcript_id, collapse = '|')))
  
  # Splitting granges by the position of enhancers and transcripts
  empP_gr_E_start <- empP_gr[start(empP_gr) < empP_gr$transcript_position]
  ep_ranges_1 <- GRanges(seqnames = rep(chr, length(empP_gr_E_start)),
                       IRanges(start = start(empP_gr_E_start),
                              end = empP_gr_E_start$transcript_position))
  
  empP_gr_T_start <- empP_gr[start(empP_gr) > empP_gr$transcript_position]
  ep_ranges_2 <- GRanges(seqnames = rep(chr, length(empP_gr_T_start)),
                         IRanges(end = start(empP_gr_T_start),
                                 start = empP_gr_T_start$transcript_position))
  
  ep_ranges <- c(ep_ranges_1, ep_ranges_2)
  # Overlapping with shuffled grb ranges, the whole range is WITHIN
  ep_ranges <- subsetByOverlaps(ep_ranges, shuffled_grb.gr, type = "within")
  ep_dist_vec <- c(ep_dist_vec, width(ep_ranges))
  ep_ranges_all <- c(ep_ranges_all, ep_ranges)
  
}

line = paste0(mean(ep_dist_vec), "\t", median(ep_dist_vec))
# write(line, file="data/random/shuffle_5085_out.txt", append=TRUE)
