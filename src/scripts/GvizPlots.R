library("stringr")
library("Gviz")
library("rtracklayer")
library("GenomicRanges")
library("biomaRt")
library("dplyr")


#setwd(dir = "02-GWAS_SNPs_inGRB")

grbs <- readRDS("data/GRBs_CNE_density_modf_kurtosis.rds")
loci_in_GRBs.df <- readRDS("data/loci_in_GRBs.df.rds")
loci_in_GRBs.gr <- readRDS("data/loci_in_GRBs.gr.rds")
snps_in_GRBs.gr <- readRDS("data/snps_in_GRBs.gr.rds")
byGRB_list <- readRDS("data/byGRB_list.rds")


# Function1
locusGviz <- function(i, dir, labels=T, byGRB_list, 
                      loci.gr, snps.gr, grbs, mart, trait) {
  
  # setting a nice pallete of Grey, Green and Purple that match
  Grey   <- "#7D949B"
  Green  <- "#9CA94E"
  Purple <- "#A74F7C" 
  
  #  BREAKPOINT: If the GRB does not have any enhancers skip it
  if(nrow(byGRB_list[[i]][[3]])==0) {
    return(print(paste("No Gviz plot for", names(byGRB_list)[i], "no enhancers found.", sep=" ")))
  }
  
  
  current_chr <- unlist(str_split(rownames(byGRB_list[[i]][[3]][1]), ":"))[1]
  
  # Saving GRB start and end pos for (later) checking if current_GRB matches byGRB
  GRB_start <- str_split(names(byGRB_list)[i], "_", simplify = TRUE)[2]
  GRB_end <- str_split(names(byGRB_list)[i], "_", simplify = TRUE)[3]
  
  s <- sapply(rownames(byGRB_list[[i]][[3]]), function(x){ as.numeric(unlist(str_split(x, "[:-]"))[2]) }) # starts of all enhancers
  e <- sapply(rownames(byGRB_list[[i]][[3]]), function(x){ as.numeric(unlist(str_split(x, "[:-]"))[3]) }) # ends of all enhancers
  range_of_enhancers <- GRanges(paste0(current_chr, ":", min(s), "-", max(e))) 
  
  current_loci <- subsetByOverlaps(loci.gr, range_of_enhancers)
  
  
  # BREAKPOINT: If the locus does not overlap any enhancers
  if(length(current_loci) == 0){
    return(print(paste("No Gviz plot for", names(byGRB_list)[i],
                       ". Trait", trait, "no loci found.", sep=" ")))
  }
  
  
  
  names(current_loci) <- current_loci$SNPS
  current_GRB <- subsetByOverlaps(grbs, current_loci)
  
  
  
  # BREAKPOINT: Checking if current_GRB contains exactly one GRB
  # and that it matches the coordinates of the byGRB_list
  if(length(current_GRB) != 1){
    current_GRB <- current_GRB[start(current_GRB) == GRB_start]
  }
  # matching starts 
  if(start(current_GRB) != as.numeric(GRB_start)){
    print("Start coordinates of GRBs do not match")
    return(message("Whoops, the current_GRB does not match byGRBlist"))
  }
  
  
  # Plot ideogram and scale
  ## Ideogram, itrack is the chromosome on top
  itrack <- IdeogramTrack(genome="hg19", chromosome=current_chr)
  
  ## GenomeAxis track, Gtrack is the track with genomic ranges
  Gtrack <- GenomeAxisTrack(name=as.character(current_chr))
  
  # GRB track 
  if(str_length(current_GRB$targets)==0) {
    grbtrack <- AnnotationTrack(current_GRB, name="GRB human:dog", 
                                id="No target",
                                fill=Grey, col=Grey,
                                genome = "hg19")
  }else{
    grbtrack <- AnnotationTrack(current_GRB, name="GRB human:dog", 
                                id=current_GRB$targets, fontcolor.feature="black", 
                                fill=Grey, col=Grey,
                                genome = "hg19")
    
    # Anjas notes - Maybe remove this?
    # really complicated way of colouring blocks into desired colour
    #feature(grbtrack) <- rep("bla", length(grbtrack))
    #displayPars(gtrack) <- list(showId=T, col=NULL, bla="red")
    
    if(labels) {
      displayPars(grbtrack) <- list(featureAnnotation="id")
    }
    # add other 6 GRB tracks - make it a grouped trackgro
  }
  
  
  # plot SNPs and GWAS-associated locus
  current_SNPs <- subsetByOverlaps(snps.gr, current_loci)
  
  
  
  
  snpTrSize <-  ifelse(length(current_SNPs)<3, 1, 
                       ifelse(length(current_SNPs)<5, 2,
                              ifelse(length(current_SNPs)<4, 3,
                                     ifelse(length(current_SNPs)<10, 4, 
                                            ifelse(length(current_SNPs)<20, 5, 6)))))
  
  
  # If current_loci has two SNPs that are identical but implicated by different 
  # studies it will appear two times. 
  # current_loci <-
  #   current_loci %>%
  #   dplyr::filter(!duplicated(SNPS))
  
  # Try to filter same LD blocks only if its too big
  if(length(current_loci) > 20){
    current_loci <- plyranges::filter(current_loci, !duplicated(LDblock))
  }
  
  ## Locus track (ltrack) LD block of a SNP
  if(length(names(current_loci))==0) {
    ltrack <- AnnotationTrack(current_loci,
                              name="GWAS locus",
                              fill=Purple, col=Purple,
                              genome = "hg19")
    #fill=Grey, col=Grey)
  }
  else {
    
    ltrack <- AnnotationTrack(current_loci,
                              name="GWAS locus", 
                              fill=Purple, col=Purple,
                              #fill=Grey, col=Grey,
                              id=current_loci$GWAS_GENES,
                              fontcolor.feature="black",
                              genome = "hg19")
    #group=names(current_loci), groupAnnotation="group") # this throws a single genome warning - ignore!
    #displayPars(ltrack) <- list(showId=T)
    if(labels) {
      displayPars(ltrack) <- list(featureAnnotation="id")
    }
  }
  
  
  lociTrSize <-  ifelse(length(current_loci)<3, 1, 
                        ifelse(length(current_loci)<5, 2,
                               ifelse(length(current_loci)<4, 3,
                                      ifelse(length(current_loci)<10, 4, 
                                             ifelse(length(current_loci)<20, 5, 6)))))
  
  
  # SNP track (strack)
  strack <- AnnotationTrack(current_SNPs, name="GWAS SNPs", 
                            stacking="full",
                            group=current_SNPs$SNPS, #groupAnnotation="group",
                            fill=Purple, col=Purple,
                            genome = "hg19")
  if(labels) {
    displayPars(strack) <- list(showId=T)
  }
  
  
  
  # Enhancer track (etrack)
  splt <- str_split(rownames(byGRB_list[[i]][[3]]), "[:-]")
  e <- GRanges(seqnames=sapply(splt, function(x){ x[1] }),
               ranges=IRanges(start=as.integer(sapply(splt, function(x) { x[2] })), 
                              end=as.integer(sapply(splt, function(x){ x[3] }))))
  ## consider only ploting enhancers that overlap loci
  e_sub <- subsetByOverlaps(e, current_loci)
  
  ## BREAKPOINT: if e_sub is empty it means the loci do not overlap any enhancers
  if(length(e_sub)==0){
    return(print("No enhancers were found in the vicinity of the SNP."))
  }
  
  
  etrack <- AnnotationTrack(e_sub, name="CAGE enhancers", 
                            stacking="squish", # default, used to be 'dense', but even in 'full' it's not expanded
                            genome = "hg19",
                            group=paste0("Enh", 1:length(e_sub)), groupAnnotation="group",
                            #fill=Grey, col=Grey) # pink is brewer.pal(8, "PiYG")[1]
                            fill=Purple, col=Purple) # pink is brewer.pal(8, "PiYG")[1]
  if(labels) {
    displayPars(etrack) <- list(showId=T) # turn on group annotation
  }
  
  # resizing track with for enhancer track
  enhTrSize <-  ifelse(length(e_sub)<4, 1, 
                       ifelse(length(e_sub)<7, 2,
                              ifelse(length(e_sub)<10, 3,
                                     ifelse(length(e_sub)<16, 4, 
                                            ifelse(length(e_sub)<25, 6, 8)))))
  
  
  # Gene track (gtrack)  
  splt <- str_split(colnames(byGRB_list[[i]][[3]]), ":")
  tr_in_locus <- sapply(splt, function(x) { return(x[2]) }) # transcripts to keep for plotting
  
  
  
  gtrack <- BiomartGeneRegionTrack(start=start(current_GRB), 
                                   end=end(current_GRB),
                                   chromosome=current_chr, 
                                   biomart=mart, 
                                   name="Genes",
                                   genome = "hg19")
  
  # only genes that are in byGRB_list (FANTOM5 transcripts)
  gtrack <- gtrack[transcript(gtrack) %in% tr_in_locus]
  
  if(length(gene(gtrack))==0){
    print(paste("error:",names(byGRB_list)[i]," Does not have genes in gtrack"))
    return(message("Whoops, something went wrong"))
  }
  
  
  # edition: added unlist here
  # Adding features to gtrack for different coloring scheme
  Gviz::feature(gtrack) <- unlist(sapply(gene(gtrack), function(g) { 
    byGRB_list[[i]][[1]]$category[byGRB_list[[i]][[1]]$ens==g]
  }))
  Gviz::feature(gtrack) <- str_replace_all(Gviz::feature(gtrack), "-", "")
  
  # Number of genes that are actually ploted (corresponds to number of symbols)
  # gene_number <- length(unique(symbol(gtrack))) 
  gene_number <- length(symbol(gtrack)) 
  
  # resizing track with for genes
  gTrSize <-  ifelse(gene_number<3, 1, 
                     ifelse(gene_number<5, 2,
                            ifelse(gene_number<8, 3,
                                   ifelse(gene_number<16, 5, 
                                          ifelse(gene_number<25, 7, 
                                                 ifelse(gene_number<37, 9, 11))))))
  
  
  # using 
  palette_primary_desktop <- c("#4E7CD9", "#153259", "#285911", "#D9A60D", "#D93425")
  col_palette <- palette_primary_desktop
  # Adding colors to genes by their category
  displayPars(gtrack) <- list(showId=T, col=NULL,
                              GWASout=col_palette[3], GWAStrgt=col_palette[5], GRBtrgt=col_palette[1], 
                              GRBbyst=col_palette[2], GWASbyst=col_palette[4]) # turn on gene symbols, no dashes allowed
  
  # CNE density track (CNEdtrack)
  #CNEdensfile <- "~/Documents/grbs/HCNE_density_hg19_mm10_96pc_50col.wig.gz"
  CNEdensfile <- "data/HCNE_density_hg19_canFam3_98pc_50col.wig.gz"
  CNEdtrack <- DataTrack(range=CNEdensfile, name="hg19_canFam3\n96/50",
                         genome="hg19", chromosome=current_chr, 
                         type="polygon", col.mountain="azure3",
                         fill.mountain=c("azure3", "azure3"))
  
  # plot directionality index
  dirfile <- "data/directionality_homer_neural.bigWig"
  dirtrack <- DataTrack(range=dirfile, name="Homer_DI_neural",
                        genome="hg19", chromosome=current_chr,
                        type="mountain", baseline=0, lty.mountain=0,
                        fill.mountain=adjustcolor(c("darkgreen", "darkred"), alpha.f = 0.7),
                        col.mountain="black") # baseline line colour
  
  dirTrSize <- 1.5 # size of the dir. index track
  # high values look noisy if not smoothed, low ones throw an error if smoothed
  # so only smooth of directionality index is above this arbitrary threshold
  DI_neural <- rtracklayer::import.bw(dirfile)
  x <- subsetByOverlaps(DI_neural, current_GRB)
  
  #adjusting span for loess to plot
  if(min(x$score)> -150 | max(x$score)<150) {
    displayPars(dirtrack)$span <- 0.25
  } # default $span=0.2
  
  if(min(x$score)< -300 & max(x$score)> 300) {
    displayPars(dirtrack)$span <- 0.15
    dirTrSize <- 2
  } # default $span=0.2
  
  #does not change much for where signal is strong, fixes span issue where low signal
  displayPars(dirtrack)$family="gaussian"
  
  ## Directionality index end
  
  #entire regions of all genes are to be covered
  xmin <- min(start(gtrack), start(current_GRB), start(current_loci))
  xmax<- max(end(gtrack), end(current_GRB), end(current_loci))
  
  # extending the view to +5% on each side of the GRB/LD/gene region
  w <- xmax-xmin
  expand <- 0.1*w/2 
  xmin <- ceiling(max(xmin-expand, 0))
  xmax <- floor(xmax+expand)
  
  # plot in a larger image if there are many enhancers/genes
  if(gTrSize==11 | enhTrSize==8) {
    wid=1000
    hei=1000
  }else{
    wid=750
    hei=750
  }
  
  
  
  
  png(filename=paste0(dir, str_replace_all(names(byGRB_list)[i], ":", "_"), "_Gviz.png"), 
      width=wid, height=hei)
  #pdf(file=paste0(dir, str_replace_all(names(byGRB_list)[i], ":", "_"), "_Gviz.pdf"),
  #    width=10, height=10)
  plotTracks(list(itrack, 
                  Gtrack,
                  grbtrack,
                  CNEdtrack,
                  dirtrack,
                  ltrack,
                  strack,
                  etrack,
                  gtrack),
             from=xmin, to=xmax, 
             sizes=c(1, #itrack 
                     1, #Gtrack
                     1, #grbtrack
                     1.5, #CNE track
                     dirTrSize,
                     lociTrSize,
                     snpTrSize,
                     enhTrSize,
                     gTrSize),
             background.title="azure3", fontcolor.title="black")
  dev.off()
}

# Creating directory for plots
gviz_dir <-"data/gviz_plots"
if(!file.exists(gviz_dir)){
  
  dir.create(gviz_dir)
}

# Adding a column to grbs that correspond to targets
grbs$targets <- grbs$geneName

# actual plotting, saved in gviz_dir
loci_in_GRBs.gr$name <- loci_in_GRBs.gr$SNPS # need to have a $name column for gviz plotting


# preparing mart to extract gene objects from
mart <- useMart(dataset="hsapiens_gene_ensembl", biomart="ENSEMBL_MART_ENSEMBL", host="https://feb2014.archive.ensembl.org")

for(trait in unique(loci_in_GRBs.df$MAPPED_TRAIT)[]){
  trait_dir <- file.path(gviz_dir,str_replace_all(trait," |/","_"))
  # Creates a dir if it does not exist
  if(!file.exists(trait_dir)){dir.create(trait_dir)}
  # to be done on each locus
  message(trait)
  # Including only loci in trait
  # Including only loci in trait
  loci_in_GRB_trait <- loci_in_GRBs.gr[loci_in_GRBs.gr$MAPPED_TRAIT == trait,]
  # Remoivng duplicated SNPS per trait/locus
  loci_in_GRB_trait <- loci_in_GRB_trait[!duplicated(loci_in_GRB_trait$SNPS),]
  
  grbs_to_loop <-
    which( names(byGRB_list) %in% str_replace_all(loci_in_GRB_trait$grb.coords, pattern = ":|-", replacement = "_"))
  
  for (i in grbs_to_loop) {
    
    print(paste("Plotting Gviz for", names(byGRB_list)[i], sep=" "))
    tryCatch(locusGviz(i, 
                       dir=paste0(trait_dir,"/",str_replace_all(trait," |/","_"),"_"),
                       labels=T, 
                       byGRB_list=byGRB_list, 
                       # Selecting only loci with the trait
                       loci.gr= loci_in_GRB_trait, 
                       # Selecting only SNPs with the trait
                       snps.gr=snps_in_GRBs.gr[snps_in_GRBs.gr$SNPS %in%
                                                 loci_in_GRB_trait$SNPS,],
                       grbs=grbs, 
                       mart=mart,
                       trait = trait),
             error=function(e){cat("ERROR :",conditionMessage(e), "\n",
                                   trait, i)})
  }
}