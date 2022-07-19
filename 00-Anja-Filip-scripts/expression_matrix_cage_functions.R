# functions, unified all for this folder (TO-DO: list for each which script uses it to filter out)
# packages and constants to be loaded
library(rtracklayer)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(corrplot)
library(plyr)
library(Gviz)
library(RColorBrewer)
library(beanplot)
options(scipen=999)

assemblies <- c("canFam3", "mm10", "monDom5", "galGal4", "lepOcu1", "danRer7", "tetNig2")

palette_primary_desktop <- c("#4E7CD9", "#153259", "#285911", "#D9A60D", "#D93425")


# will go through all files in dir with pattern (assumes these are .bed files)
# every range in loci.gr will get TRUE if it overlaps any of the ranges in any of the .bed files in dir
# outputs a logical vector, one value for each range in loci.gr
is_overlapping_any_GRB <- function(loci.gr, dir, file_pattern="processed.threshold.grbs.bed") {

  v <- sapply(list.files(dir, pattern=file_pattern), function(f) {
    
      gr <- import.bed(paste0(dir, f))
      overlapsAny(loci.gr, gr)
  })
  #return(apply(v, 1, any)) # if TRUE for any of the pairwise comparisons, will return TRUE
}


# will go through all files in dir with pattern (assumes these are .bed files)
# every range in loci.gr will get TRUE if it overlaps any of the ranges in any of the .bed files in dir
# outputs a matrix, one row for each range, one column per file in dir
is_overlapping_GRBs_all <- function(query_gr, dir, file_pattern="processed.threshold.grbs.bed") {
  
  v <- sapply(list.files(dir, pattern=file_pattern), function(f) {
    
    gr <- import.bed(paste0(dir, f))
    overlapsAny(query_gr, gr)
  })
  colnames(v) <- paste0("in_grb_", assemblies)
  return(v)
}


#which.id can be "ens" or "hgnc"
extract_GRBtargets_Ge_all <- function(loci.gr, which.id="ens", name.sep=":") {
  for (a in assemblies) {
    print(paste0("[make_loci_gr] Extracting target genes for GRBs human vs. ", a))
    
    grb <- import.bed(grep(a, list.files("~/Documents/grbs/", pattern="processed.threshold", full.names=T), value=T))
    
    # gr is GRanges of target genes predicted by Ge's ML method, lifted over to hg19
    if(a=="mm10") {
      t <- read.table(paste("~/Documents/grbs/threshold_for_target_reanalysis/new_", a, ".tsv", sep=""), skip=1, stringsAsFactors=F)
    }
    else{
      t <- read.table(paste("~/Documents/grbs/target_genes_Ge/hg19.", a, ".tsv", sep=""), skip=1, stringsAsFactors=F)
    }
    splt <- matrix(unlist(str_split(t[,6], pattern="[:-]")), ncol=3, byrow=T)
    gr <- GRanges(seqnames=splt[,1],
                  ranges=IRanges(start=as.numeric(splt[,2]),
                                 end=as.numeric(splt[,3])),
                  gene_symb=t[,2],
                  gene_ens=t[,1]
                  )
    # t[,7] is Anja's 'p>0.5 no paralogue condition' prediction, 
    # t[,5] is the original Ge's normalisation
    if( a!="mm10") {
      gr$prediction <- t[,5]
    }
    else {
      gr$prediction <- t[,7]
    }
    gr <- gr[gr$prediction=="target"]
    
    ovlp1 <- findOverlaps(loci.gr, grb) # first identify which grb locus begins to 
    
    targets <- character(length(loci.gr))
    for(i in 1:length(ovlp1)) {
      
      ovlp2 <- findOverlaps(grb[subjectHits(ovlp1)[i]], gr) # this yields list of target genes for that GRB
      
      if(which.id=="ens") {
        targets[queryHits(ovlp1)[i]] <- paste(gr$gene_ens[subjectHits(ovlp2)], collapse=",") 
        #comma-separate them and assign to right locus
      }
      else if(which.id=="hgnc") {
        targets[queryHits(ovlp1)[i]] <- paste(gr$gene_symb[subjectHits(ovlp2)], collapse=",") 
        #comma-separate them and assign to right locus
      }
      else {
        print("[extract_GRBtargets_Ge_all] which.id needs to be 'ens' or 'hgnc'!")
        return(NA)
      }
    }
    
    loci.gr$targets <- targets
    colnames(values(loci.gr))[ncol(values(loci.gr))] <- paste("targets_ge_", which.id, "_", a, sep="")
  }
  loci.gr$name <- separate_non_empty_per_row(as.data.frame(mcols(loci.gr)), ":", which.id)
  colnames(values(loci.gr))[ncol(values(loci.gr))] <- paste("GRBname_", which.id, sep="")
  return(loci.gr)
}


prepare_hg19_genes_Ensembl_from_UCSC <- function() {
  library(GenomicFeatures)
  txdbEnsembl <- makeTxDbFromUCSC("hg19", "ensGene")
  # consider doing biomaRt extraction from archives

  genes <- genes(txdbEnsembl)
  
  library(BSgenome.Hsapiens.UCSC.hg19)
  seqlengths(genes) <- seqlengths(Hsapiens)[names(seqlengths(genes))]
  return(unique(genes))
}


prepare_hg19_genes_Ensembl_from_biomaRt <- function() {
  library(biomaRt)
  library(BSgenome.Hsapiens.UCSC.hg19)
  
  bm.annotations = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", 
                                      "start_position", "end_position", "strand"), 
                         mart=human, uniqueRows=TRUE)
  
  bm.annotations <- bm.annotations[which(str_length(bm.annotations$chromosome_name)<3),]
  genes <- GRanges(seqnames=bm.annotations$chromosome_name,
                   ranges=IRanges(start=bm.annotations$start_position, end=bm.annotations$end_position),
                   strand=ifelse(bm.annotations$strand==1, "+", ifelse(bm.annotations$strand==-1, "-", "*")),
                   hgnc_symbol=bm.annotations$hgnc_symbol)
  seqlevelsStyle(genes) <- "UCSC" # bm.annotations is 'NCBI'
  seqlengths(genes) <- seqlengths(Hsapiens)[names(seqlengths(genes))]
  
  names(genes) <- bm.annotations$ensembl_gene_id
  return(unique(genes))
}


load_all_human_FANTOM5_samples <- function() {
  
  #library(CAGEr)
  # load all samples
  #f <- "./FANTOM5/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt"
  # first column as row names, then 2-7 are protID (annotation), then tpm starts
  
  #cmd <- paste0("grep -v '##' ", f, " | head -n 1 | awk '{print NF}'")
  #ncolumns <- as.integer(system(cmd, intern=T)) # to define number of columns
  # there are 1837 rows with comments, then colnames row, then '01STAT:MAPPED' then '02STAT:NORM_FACTOR' line
  # first 7 columns (rownames+6) are identifiers which protein peak maps to, then TPMs start
  #annot <- read.table(file=f, sep="\t", row.names=1, nrows=202000, skip=1840, stringsAsFactors=F)
  #namesOfColumns <- read.table(file=f, sep="\t", row.names=1, nrows=1, skip=1837, stringsAsFactors=F)
  
  # first 6 cols are annotation, no need for this at this stage
  #annot.6cols <- annot[,1:6]
  
  #tpm.df <- annot[,7:ncol(annot)]
  #colnames(tpm.df) <- unname(namesOfColumns[7:length(namesOfColumns)])
  #saveRDS(tpm.df, file="FANTOM5_tpm_annot.df.rds")
  #rm(annot)
  tpm.df <- readRDS("FANTOM5_tpm_annot.df.rds")
  
  # make GRanges out of it
  a <- unlist(str_split(rownames(tpm.df), "[:,..]")) 
  m <- matrix(a, ncol=5, byrow=T) # third column are empty strings because of '..'
  gr <- GRanges(seqnames=m[,1],
                ranges=IRanges(start=as.integer(m[,2]), 
                               end=as.integer(m[,4])),
                strand=m[,5])
  mcols(gr) <- tpm.df
  return(gr)
}


load_all_human_FANTOM5_samples_enhancer_expression <- function() {
  library(CAGEr)
  # load all samples
  #f <- "~/Documents/SchizophreniaSNPs/FANTOM5/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt"
  
  # there is a line with sample names and 65423 lines with per-enhancer TPM-normalised expression
  # first column is peak ID then expression values
  #annot <- read.table(file=f, sep="\t", row.names=1, nrows=66000, skip=1, stringsAsFactors=F)
  #namesOfColumns <- read.table(file=f, sep="\t", row.names=1, nrows=1, stringsAsFactors=F)
  #colnames(annot) <- unname(namesOfColumns)

  #saveRDS(annot, file="FANTOM5_tpm_annot_enhancers.df.rds")
  tpm.df <- readRDS("FANTOM5_tpm_annot_enhancers.df.rds")
  
  # make GRanges out of it
  a <- unlist(str_split(rownames(tpm.df), "[:-]")) 
  m <- matrix(a, ncol=3, byrow=T) # third column are empty strings because of '..'
  gr <- GRanges(seqnames=m[,1],
                ranges=IRanges(start=as.integer(m[,2]), 
                               end=as.integer(m[,3])))
  mcols(gr) <- tpm.df
  return(gr)
}

load_merged_human_FANTOM5_samples_enhancer_expression <- function() {
  # library(CAGEr)
  # # load all samples
  # f <- "~/Documents/snpsInGrbs/precalcData/FANTOM5/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt"
  # 
  # # there is a line with sample names and 65423 lines with per-enhancer TPM-normalised expression
  # # first column is peak ID then expression values
  # annot <- read.table(file=f, sep="\t", row.names=1, nrows=66000, skip=1, stringsAsFactors=F)
  # namesOfColumns <- read.table(file=f, sep="\t", row.names=1, nrows=1, stringsAsFactors=F)
  # colnames(annot) <- unname(namesOfColumns)
  # 
  # # now do the same sample pooling as in promoters
  # # this also eliminates all extra samples, only tissues and primary cells are included in this list
  # biolRepList <- readRDS("biolRepList.rds")
  # logannot <- log_transform(annot, 10)
  # logannot_merged <- merge_replicate_log_tpms(biolRepList, logannot) # this step requires logged data
  # annot_merged <- rev_log_transform(logannot_merged, 10) # go back from logs
  #saveRDS(annot_merged, file="FANTOM5_merged_tpm_annot_enhancers.df.rds")
  
  tpm.df <- readRDS("FANTOM5_merged_tpm_annot_enhancers.df.rds")
  
  # make GRanges out of it
  a <- str_split(rownames(tpm.df), "[:-]", simplify=T)
  gr <- GRanges(seqnames=a[,1],
                ranges=IRanges(start=as.integer(a[,2]),
                               end=as.integer(a[,3])))
  mcols(gr) <- tpm.df
  return(gr)
}



# orders mcols by colnames_ordered, and eliminates any missing from colnames_ordered
# returns orig. GRanges objects with alteres mcols
reorder_prom_colnames_by_enh <- function(gr, colnames_ordered) {
  
  # prom colnames have full names with tissue info etc., needs trimming
  # splits by dot, then extracts the sampleID one
  l <- lapply(str_split(colnames(mcols(gr)), "\\."), function(x) {
    grep("CNhs", x, value=T)
  })
  colnames(mcols(gr)) <- unlist(l)
  
  mcols(gr) <- mcols(gr)[,colnames_ordered] # will removed samples not existing in colnames_ordered, and reorder the rest
  return(gr)
}


load_all_human_FANTOM5_sample_names <- function() {
  f <- "~/Documents/snpsInGrbs/precalcData/FANTOM5/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt"
  namesOfColumns <- read.table(file=f, sep="\t", row.names=1, nrows=1, skip=1837, stringsAsFactors=F)
}


# requires a GRanges with promoters and another one with CAGE peaks, where mcols are one sample per column
# out puts a data.frame with a row per promoter and a column per FANTOM5 sample with summed tpms per sample
old_sum_peak_signal_per_promoter <- function(promoters.gr, fantom5.tpm.gr) {
  
  ovlp <- findOverlaps(promoters.gr, fantom5.tpm.gr)
  
  signal_per_promoter <- sapply(unique(queryHits(ovlp)), function(i) {
    #i is index of which ps entry we are working on
    which_subject <- subjectHits(ovlp[queryHits(ovlp)==i])
    per_sample_summed_signal <- apply(as.data.frame(mcols(fantom5.tpm.gr[which_subject])), 2, sum)
  })
  mcols_to_fill <- t(signal_per_promoter)
  rownames(mcols_to_fill) <- unique(queryHits(ovlp))
  
  df <- matrix(rep(0, length(promoters.gr)*ncol(mcols(fantom5.tpm.gr))), 
               nrow=length(promoters.gr), ncol=ncol(mcols(fantom5.tpm.gr)))
  rownames(df) <- names(promoters.gr)
  
  for(i in unique(queryHits(ovlp))) {
    if(i %in% rownames(mcols_to_fill)) {
      df[i,] <- mcols_to_fill[rownames(mcols_to_fill)==i,]
    }
  }
  return(df)
}


# requires a GRanges with promoters and another one with CAGE peaks, where mcols are one sample per column
# out puts a data.frame with a row per promoter and a column per FANTOM5 sample with summed tpms per sample
# a gene can have multiple promoters, thus rownames are ENSG:ENST format
sum_peak_signal_per_promoter <- function(promoters.gr, fantom5.tpm.gr) {
  
  ovlp <- findOverlaps(promoters.gr, fantom5.tpm.gr)
  
  signal_per_promoter <- sapply(unique(queryHits(ovlp)), function(i) {
    #i is index of which pr entry we are working on
    if (i%%100==0) {
      print(i)
    }
    which_subject <- subjectHits(ovlp[queryHits(ovlp)==i])
    per_sample_summed_signal <- apply(as.data.frame(mcols(fantom5.tpm.gr[which_subject])), 2, sum)
  })
  mcols_to_fill <- t(signal_per_promoter)
  
  df <- matrix(rep(0, length(promoters.gr)*ncol(mcols(fantom5.tpm.gr))), 
               nrow=length(promoters.gr), ncol=ncol(mcols(fantom5.tpm.gr)))
  df[unique(queryHits(ovlp)),] <- mcols_to_fill
  rownames(df) <- paste(promoters.gr$ensembl_gene_id, names(promoters.gr), sep=":")
  colnames(df) <- colnames(mcols_to_fill)

  return(df)
}


# input is a data frame, output a vector of nrow length
# collates all non-empty strings using colon as a separator
separate_non_empty_per_row <- function(df, sep=":", pattern="") {
  v <- apply(df[,grep(pattern, colnames(df), value=T)], 1, function(x) {
    x <- unlist(str_split(x, ","))
    x <- x[which(str_length(x)>0)]
    y <- paste(unique(x), collapse=sep)
    return(str_replace_all(y, paste0(sep, sep), sep)) # to remove "::"
  })
  return(v) 
}


#input a vector of gwas-associated genes
gwas_gene_cleanup <- function(gwas, sep=":") {
  gwas <- str_replace(gwas, " \\*", "")
  gwas <- str_replace(gwas, " \\(micro-RNA\\)", "")
  gwas <- str_replace(gwas, "Locus too broad", "")  
  gwas <- str_replace_all(gwas, " ", sep)
}


# current_GRB is a GRange of current GRB to be processed
# current_loci needs to have a GRanges with ranges of the loci overlapping that GRB
# and (colon delimited) string of GWAS associated in $gwas_genes field, and (colon delimited) string of target genes in $GRBname field (all target s predicted in all 7 pairwise comparisons)
# outputs a df with 3 columns, one row for each gene in that GRB, 
# plus any GWAS-associated genes with loci overlapping that GRB that fall outside of the GRB
# first column is EnsemblID, secon id hgnc, third can be one of 5 classes
# (i) 'GWAS-out' for a GWAS-associated gene falling outside of GRB boundaries
# (ii) 'GWAS-trgt' for GWAS-associated gene that is a predicted target in that GRB
# (iii) 'GWAS-byst' for all other GWAS-associated genes
# (iv) 'GRB-trgt' for all other target genes predicted in that GRB (not previously associated)
# (v) 'GRB-byst' all other genes in that GRB (not previously associated and not target genes)
per_GRB_genes_classification <- function(current_GRB, loci.gr) {
  
  current_loci <- subsetByOverlaps(loci.gr, current_GRB)
  current_chr <- as.character(seqnames(current_loci))[1]
  g <- subsetByOverlaps(genes.hg19.ensembl, current_GRB)
  
  # add GRB genes
  genes.df <- data.frame(ens=names(g),
                         strt=start(g),
                         symbol=g$hgnc_symbol,
                         #category=character(length(g)), # category is added in GWAS step
                         stringsAsFactors=F)
  rownames(genes.df) <- as.character(g)
  
  # add GWAS genes
  if(any(str_length(current_loci$gwas_genes)>0)) {
    symb <- unique(unlist(str_split(current_loci$gwas_genes, ",")))
    symb <- symb[-which(symb=="")]
    gwas <- data.frame(symbol=symb)
    
    gwas$ensembl <- sapply(gwas$symbol, function(x) {
      # it can be one-to-one, one-to-many and one-to-none mapping
      a <- names(genes.hg19.ensembl[genes.hg19.ensembl$hgnc_symbol==x & as.character(seqnames(genes.hg19.ensembl))==current_chr])
      if(length(a)==1) { 
        return(a)
      }
      else if (length(a)<1) {
        return(NA)
      }
      else {
        
        print("Oooops! Found multiple ensIDs for this gene symbol:")
        print(a)
        return(NA)
      }
    })
    gwas <- na.omit(gwas)
    gwas$in_GRB <- ifelse(gwas$ensembl %in% genes.df$ens, "GWAS-in", "GWAS-out")
    gwas$start <- sapply(gwas$symbol, function(h) {
      start(genes.hg19.ensembl[genes.hg19.ensembl$hgnc_symbol==h])
    })
    
    genes.df <- merge(genes.df, gwas, by.x=c(1,3,2), by.y=c(2,1,4), all=T)
    colnames(genes.df)[4] <- "category"
    genes.df$category[which(is.na(genes.df$category))] <- ""
  }
  else {
    genes.df$category <- ""
  }
  
  # add target genes
  #target <- data.frame(ensembl=unlist(str_split(current_loci$targets_ge_ens_mm10, ",")), stringsAsFactors=F)
  target <- data.frame(symbol=unlist(str_split(current_loci$GRBtargets, ",")), 
                       stringsAsFactors=F)
  target$category <- sapply(target$symbol, function(x) {
    cat <- ifelse(genes.df$category[genes.df$symbol==x]=="GWAS-in", "GWAS-trgt", "GRB-trgt")[1]
  })
  target <- unique(target)
  
  for (i in 1:nrow(target)) {
    genes.df$category[genes.df$symbol==target$symbol[i]] <- target$category[i]
  }
  
  # all others are bystanders
  genes.df$category <- ifelse(genes.df$category=="", 
                              "GRB-byst",
                              ifelse(genes.df$category=="GWAS-in", 
                                     "GWAS-byst",
                                     genes.df$category))
  genes.df <- genes.df[order(genes.df$strt),]
  genes.df$strt <- NULL
  return(na.omit(unique(genes.df)))
}


# loads a per-chromosome rho file and returns a GRanges with all non-NA rho values
# there will be one entry per each enhancer-gene pair, range will be for enhancer
load_per_chr_rho_file <- function(current_chr) {
  
  f <- paste0("C:/Users/Anja/Documents/Projects/SchizophreniaSNPs/rho_files/cageTSSs_without_zeroes/rho_fixed_CAGEenhPermissive_hg19Ens_6.4Mb_", as.character(current_chr), ".txt")
  print(paste("[load_per_chr_rho_file] Loading", f, sep=" "))
  rho_values <- read.table(f, sep="\t", stringsAsFactors=F, skip=1)
  
  # remove the rows with no rho value, it should be penultimate column
  toremove <- which(is.na(rho_values[,ncol(rho_values)-1]))
  rho_values <- rho_values[-toremove,]
  
  colnames(rho_values) <- scan(f, what=character(), sep="\t", nlines=1)
  
  splt <- matrix(unlist(str_split(rho_values$enh, pattern="[:-]")), ncol=3, byrow=T)
  rho_gr <- GRanges(seqnames=rho_values$chr,
                    ranges=IRanges(start=as.integer(splt[,2]), end=as.integer(splt[,3])),
                    gene=rho_values$gene,
                    rho=rho_values$rho,
                    pval=rho_values$pval)
}


load_per_chr_mi_file <- function(current_chr) {
  
  f <- paste0("C:/Users/Anja/Documents/Projects/SchizophreniaSNPs/mi_files/mi_CAGEenhPermissive_hg19Ens_6.4Mb_", as.character(current_chr), ".txt")
  print(paste("[load_per_chr_mi_file] Loading", f, sep=" "))
  mi_values <- read.table(f, sep="\t", stringsAsFactors=F, skip=1)
  mi_values <- na.omit(mi_values)
  colnames(mi_values) <- scan(f, what=character(), sep="\t", nlines=1)
  
  splt <- matrix(unlist(str_split(mi_values$enh, pattern="[:-]")), ncol=3, byrow=T)
  mi_gr <- GRanges(seqnames=mi_values$chr,
                    ranges=IRanges(start=as.integer(splt[,2]), end=as.integer(splt[,3])),
                    gene=mi_values$gene,
                    rho=mi_values$rho)
}


load_per_chr_empP_file <- function(current_chr) {
  
  f <- paste0("C:/Users/Anja/Documents/Projects/snpsInGrbs/precalcData/empP_all_files/cageTSSs_pooledFANTOMsamples/empP_CAGEenhPermissive_hg19Ens_6.4Mb_", as.character(current_chr), ".txt")
  print(paste("[load_per_chr_empP_file] Loading", f, sep=" "))
  empP_values <- read.table(f, sep="\t", stringsAsFactors=F, skip=1)
  #empP_values <- na.omit(empP_values)
  colnames(empP_values) <- scan(f, what=character(), sep="\t", nlines=1)
  
  splt <- str_split(empP_values$enh, pattern="[:-]", simplify=T)
  empP_gr <- GRanges(seqnames=empP_values$chr,
                     ranges=IRanges(start=as.integer(splt[,2]), end=as.integer(splt[,3])),
                     gene=empP_values$gene,
                     empP=empP_values$empP)
}


per_GRB_correlation_matrix <- function(rho_gr, current_GRB, genes_to_include) {
  
  library(reshape2)
  current_loci <- subsetByOverlaps(loci.gr, current_GRB)
  enhancers_in_loci <- subsetByOverlaps(rho_gr, current_loci) # looking at correlations of (permissive) enhancers falling within LD loci only
  # enhancers_in_loci <- enhancers_in_loci[enhancers_in_loci$gene %in% genes_to_include] 
  # only keep the scores to relevant genes
  enhancers_in_loci <- enhancers_in_loci[grep(paste(genes_to_include, collapse="|"),
                                              enhancers_in_loci$gene)]
  
  if(length(enhancers_in_loci)==0) {
    return(data.frame())
  }
  else {
    df <- data.frame(as.character(enhancers_in_loci), enhancers_in_loci$gene, enhancers_in_loci$rho, enhancers_in_loci$pval)
    df_tab <- dcast(df[,1:3], as.character.enhancers_in_loci. ~ enhancers_in_loci.gene, value.var="enhancers_in_loci.rho")
    rownames(df_tab) <- df_tab[,1]
    df_tab[,1] <- NULL
    return(df_tab)  
  }
}


per_GRB_mi_matrix <- function(mi_gr, current_GRB, genes_to_include) {
  
  library(reshape2)
  current_loci <- subsetByOverlaps(loci.gr, current_GRB)
  enhancers_in_loci <- subsetByOverlaps(mi_gr, current_loci) # looking at correlations of (permissive) enhancers falling within LD loci only
  # enhancers_in_loci <- enhancers_in_loci[enhancers_in_loci$gene %in% genes_to_include] 
  # only keep the scores to relevant genes
  enhancers_in_loci <- enhancers_in_loci[grep(paste(genes_to_include, collapse="|"),
                                              enhancers_in_loci$gene)]
  
  if(length(enhancers_in_loci)==0) {
    return(data.frame())
  }
  else {
    df <- data.frame(as.character(enhancers_in_loci), 
                     enhancers_in_loci$gene, 
                     enhancers_in_loci$rho)
    df_tab <- dcast(df[,1:3], 
                    as.character.enhancers_in_loci. ~ enhancers_in_loci.gene, 
                    value.var="enhancers_in_loci.rho")
    rownames(df_tab) <- df_tab[,1]
    df_tab[,1] <- NULL
    return(df_tab)  
  }
}


per_GRB_empP_matrix <- function(empP_gr, current_GRB, genes_to_include) {
  
  library(reshape2)
  current_loci <- subsetByOverlaps(loci.gr, current_GRB)
  enhancers_in_loci <- subsetByOverlaps(empP_gr, current_loci) # looking at correlations of (permissive) enhancers falling within LD loci only
  # enhancers_in_loci <- enhancers_in_loci[enhancers_in_loci$gene %in% genes_to_include] 
  # only keep the scores to relevant genes
  enhancers_in_loci <- enhancers_in_loci[grep(paste(genes_to_include, collapse="|"),
                                              enhancers_in_loci$gene)]
  
  if(length(enhancers_in_loci)==0) {
    return(data.frame())
  }
  else {
    df <- data.frame(as.character(enhancers_in_loci), enhancers_in_loci$gene, enhancers_in_loci$empP)
    df_tab <- dcast(df[,1:3], as.character.enhancers_in_loci. ~ enhancers_in_loci.gene, value.var="enhancers_in_loci.empP")
    rownames(df_tab) <- df_tab[,1]
    df_tab[,1] <- NULL
    return(df_tab)  
  }
}


hack_loci_pc_genes <- function(pc_genes) {
  
  pc_genes <- str_replace_all(pc_genes, "CNNM2CYP17A1", "CNNM2 CYP17A1")
  pc_genes <- str_replace_all(pc_genes, "PCGF6PDCD11", "PCGF6 PDCD11")
  pc_genes <- str_replace_all(pc_genes, "TRIM8USMG5", "TRIM8 USMG5")
  pc_genes <- str_replace_all(pc_genes, "C12orf65CDK2AP1", "C12orf65 CDK2AP1")
  pc_genes <- str_replace_all(pc_genes, "MPHOSPH9OGFOD2", "MPHOSPH9 OGFOD2")
  pc_genes <- str_replace_all(pc_genes, "RILPL2SBNO1", "RILPL2 SBNO1")
  pc_genes <- str_replace_all(pc_genes, "C2orf69TYW5", "C2orf69 TYW5")
  pc_genes <- str_replace_all(pc_genes, "BAG5CKB", "BAG5 CKB")
  pc_genes <- str_replace_all(pc_genes, "CHRNA3CHRNA5", "CHRNA3 CHRNA5")
  pc_genes <- str_replace_all(pc_genes, "IREB2PSMA4", "IREB2 PSMA4")
  pc_genes <- str_replace_all(pc_genes, "NDUFA4L2NXPH4", "NDUFA4L2 NXPH4")
  pc_genes <- str_replace_all(pc_genes, "SHMT2STAC3", "SHMT2 STAC3")
  pc_genes <- str_replace_all(pc_genes, "TAC3TMEM194A", "TAC3 TMEM194A")
  pc_genes <- str_replace_all(pc_genes, "GIGYF2KCNJ13", "GIGYF2 KCNJ13")
  pc_genes <- str_replace_all(pc_genes, "ATG13CHRM4", "ATG13 CHRM4")
  pc_genes <- str_replace_all(pc_genes, "CREB3L1DGKZ", "CREB3L1 DGKZ")
  pc_genes <- str_replace_all(pc_genes, "GOLGA6L4ZSCAN2", "GOLGA6L4 ZSCAN2")
  pc_genes <- str_replace_all(pc_genes, "PLCL1RFTN2", "PLCL1 RFTN2")
  pc_genes <- str_replace_all(pc_genes, "L3MBTL2RANGAP1", "L3MBTL2 RANGAP1")
  pc_genes <- str_replace_all(pc_genes, "GLT8D1 GNL3 ITIH1 ITIH3", "GLT8D1 GNL3 ITIH1 ITIH3 ITIH4 MUSTN1 NEK4 NISCH NT5DC2 PBRM1 SMIM4 SPCS1 STAB1 TMEM110 TMEM110-MUSTN1") # missed genes while copy-pasting from .pdf
  pc_genes <- str_replace_all(pc_genes, "C16orf92DOC2A", "C16orf92 DOC2A")
  pc_genes <- str_replace_all(pc_genes, "GDPD3HIRIP3", "GDPD3 HIRIP3")
  pc_genes <- str_replace_all(pc_genes, "KCTD13MAPK3", "KCTD13 MAPK3")
  pc_genes <- str_replace_all(pc_genes, "TAOK2TBX6", "TAOK2 TBX6")
  pc_genes <- str_replace_all(pc_genes, "HAPLN4MAU2", "HAPLN4 MAU2")
  pc_genes <- str_replace_all(pc_genes, "PBX4SUGP1", "PBX4 SUGP1")
  pc_genes <- str_replace_all(pc_genes, "C1orf51C1orf54", "C1orf51 C1orf54")
  pc_genes <- str_replace_all(pc_genes, "SHISA8SMDT1", "SHISA8 SMDT1")
  pc_genes <- str_replace_all(pc_genes, "TCF20TNFRSF13C", "TCF20 TNFRSF13C")
  pc_genes <- str_replace_all(pc_genes, "CLP1CTNND1", "CLP1 CTNND1")
  pc_genes <- str_replace_all(pc_genes, "SERPING1TMX2", "SERPING1 TMX2")
  pc_genes <- str_replace_all(pc_genes, "ETF1FAM53C", "ETF1 FAM53C")
  pc_genes <- str_replace_all(pc_genes, "HSPA9KDM3B", "HSPA9 KDM3B")
  pc_genes <- str_replace_all(pc_genes, "GFOD2LCAT", "GFOD2 LCAT")
  pc_genes <- str_replace_all(pc_genes, "NUTF2PARD6A", "NUTF2 PARD6A")
  pc_genes <- str_replace_all(pc_genes, "PSKH1PSMB10", "PSKH1 PSMB10")
  pc_genes <- str_replace_all(pc_genes, "SLC12A4SLC7A6", "SLC12A4 SLC7A6")
  pc_genes <- str_replace_all(pc_genes, "THAP11TSNAXIP1", "THAP11 TSNAXIP1")
  pc_genes <- str_replace_all(pc_genes, "LRRC48MYO15A", "LRRC48 MYO15A")
  pc_genes <- str_replace_all(pc_genes, "SREBF1TOM1L2", "SREBF1 TOM1L2")
  pc_genes <- str_replace_all(pc_genes, "RCN3RRAS", "RCN3 RRAS")
  pc_genes <- str_replace_all(pc_genes, "PCDHA1PCDHA10", "PCDHA1 PCDHA10")
  pc_genes <- str_replace_all(pc_genes, "PCDHA3PCDHA4", "PCDHA3 PCDHA4")
  pc_genes <- str_replace_all(pc_genes, "PCDHA6PCDHA7", "PCDHA6 PCDHA7")
  pc_genes <- str_replace_all(pc_genes, "PCDHA9TMCO6", "PCDHA9 TMCO6")
}


# need to remove all time series (and other dodgy) samples
# tried loading from CAGEr package, but type there is "cell line", "primary cell" and "tissue"
# decided to compare against our list of .bed files downloaded in *** where time series is separate
fantom_sample_cleanup <- function() {
  
  #library(CAGEr)
  #data(FANTOM5humanSamples)
  # one way would be to only include datasets from CAGEr that have type 'cell line', 'primary cell' or tissue'
  # alternatively, all files are downloaded to mammoth at: /mnt/biggles/csc_projects/piotr/FANTOM5/Data/PrimaryData/UPDATE_023
  # obtained list of library IDs for time series ones there and keep the rest
  
  
  sample_names <- load_all_human_FANTOM5_sample_names()
#  s <- sapply(sample_names, function(x) {
#    unlist(str_split(x, "\\."))[2]
#  })
#  s <- unname(s)
  
  c <- sapply(sample_names, function(x) {
    unlist(str_split(x, "\\."))[3]
  })
  c <- unname(c)
  
  time_series <- scan("~/Documents/snpsInGrbs/precalcData/FANTOM5/FANTOM5.human.timecourse.files.txt", what=character())
  t <- grep("CNhs", unlist(str_split(time_series, "\\.")), value=T)
  return (c[!(c %in% t)])
}


# corrplot::corrplot adapted to accept 2 tl.col vectors and colours x and y labels accordingly
# if tl.pos='d' then just uses tl.col.x
corrplot_2col_labels <- function(corr, method = c("circle", "square", "ellipse", "number", 
                                                   "shade", "color", "pie"), type = c("full", "lower", "upper"), 
                                  add = FALSE, col = NULL, bg = "white", title = "", is.corr = TRUE, 
                                  diag = TRUE, outline = FALSE, mar = c(0, 0, 0, 0), addgrid.col = NULL, 
                                  addCoef.col = NULL, addCoefasPercent = FALSE, order = c("original", "AOE", "FPC", "hclust", "alphabet"), hclust.method = c("complete", "ward", "single", "average", "mcquitty", "median", "centroid"), 
                                  addrect = NULL, rect.col = "black", rect.lwd = 2, tl.pos = NULL, 
                                  tl.cex = 1, tl.col.x = "red", tl.col.y = "red", tl.offset = 0.4, tl.srt = 90, 
                                  cl.pos = NULL, cl.lim = NULL, cl.length = NULL, cl.cex = 0.8, 
                                  cl.ratio = 0.15, cl.align.text = "c", cl.offset = 0.5, addshade = c("negative", "positive", "all"), shade.lwd = 1, shade.col = "white", 
                                  p.mat = NULL, sig.level = 0.05, insig = c("pch", "p-value", 
                                                                            "blank", "n"), pch = 4, pch.col = "black", pch.cex = 3, 
                                  plotCI = c("n", "square", "circle", "rect"), lowCI.mat = NULL, 
                                  uppCI.mat = NULL, ...) 
{
  method <- match.arg(method)
  type <- match.arg(type)
  order <- match.arg(order)
  hclust.method <- match.arg(hclust.method)
  plotCI <- match.arg(plotCI)
  insig <- match.arg(insig)
  if (!is.matrix(corr) & !is.data.frame(corr)) 
    stop("Need a matrix or data frame!")
  if (is.null(addgrid.col)) {
    addgrid.col <- ifelse(method == "color" | method == "shade", 
                          "white", "grey")
  }
  if (any(corr < cl.lim[1]) | any(corr > cl.lim[2])) 
    stop("color limits should cover matrix")
  if (is.null(cl.lim)) {
    if (is.corr) 
      cl.lim <- c(-1, 1)
    if (!is.corr) 
      cl.lim <- c(min(corr), max(corr))
  }
  intercept <- 0
  zoom <- 1
  if (!is.corr) {
    if (max(corr) * min(corr) < 0) {
      intercept <- 0
      zoom <- 1/max(abs(cl.lim))
    }
    if (min(corr) >= 0) {
      intercept <- -cl.lim[1]
      zoom <- 1/(diff(cl.lim))
    }
    if (max(corr) <= 0) {
      intercept <- -cl.lim[2]
      zoom <- 1/(diff(cl.lim))
    }
    corr <- (intercept + corr) * zoom
  }
  cl.lim2 <- (intercept + cl.lim) * zoom
  int <- intercept * zoom
  if (min(corr) < -1 - .Machine$double.eps || max(corr) > 1 + 
      .Machine$double.eps) {
    stop("The matrix is not in [-1, 1]!")
  }
  if (is.null(col)) {
    col <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", 
                              "#F4A582", "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE", 
                              "#4393C3", "#2166AC", "#053061"))(200)
  }
  n <- nrow(corr)
  m <- ncol(corr)
  min.nm <- min(n, m)
  ord <- 1:min.nm
  if (!order == "original") {
    ord <- corrMatOrder(corr, order = order, hclust.method = hclust.method)
    corr <- corr[ord, ord]
  }
  if (is.null(rownames(corr))) 
    rownames(corr) <- 1:n
  if (is.null(colnames(corr))) 
    colnames(corr) <- 1:m
  getPos.Dat <- function(mat) {
    x <- matrix(1:n * m, n, m)
    tmp <- mat
    if (type == "upper") 
      tmp[row(x) > col(x)] <- Inf
    if (type == "lower") 
      tmp[row(x) < col(x)] <- Inf
    if (type == "full") 
      tmp <- tmp
    if (!diag) 
      diag(tmp) <- Inf
    Dat <- tmp[is.finite(tmp)]
    ind <- which(is.finite(tmp), arr.ind = TRUE)
    Pos <- ind
    Pos[, 1] <- ind[, 2]
    Pos[, 2] <- -ind[, 1] + 1 + n
    return(list(Pos, Dat))
  }
  Pos <- getPos.Dat(corr)[[1]]
  n2 <- max(Pos[, 2])
  n1 <- min(Pos[, 2])
  nn <- n2 - n1
  newrownames <- as.character(rownames(corr)[(n + 1 - n2):(n + 
                                                             1 - n1)])
  m2 <- max(Pos[, 1])
  m1 <- min(Pos[, 1])
  mm <- m2 - m1
  newcolnames <- as.character(colnames(corr)[m1:m2])
  DAT <- getPos.Dat(corr)[[2]]
  len.DAT <- length(DAT)
  assign.color <- function(DAT) {
    newcorr <- (DAT + 1)/2
    newcorr[newcorr == 1] <- 1 - 1e-10
    col.fill <- col[floor(newcorr * length(col)) + 1]
  }
  col.fill <- assign.color(DAT)
  isFALSE = function(x) identical(x, FALSE)
  isTRUE = function(x) identical(x, TRUE)
  if (isFALSE(tl.pos)) {
    tl.pos <- "n"
  }
  if (is.null(tl.pos) | isTRUE(tl.pos)) {
    if (type == "full") 
      tl.pos <- "lt"
    if (type == "lower") 
      tl.pos <- "ld"
    if (type == "upper") 
      tl.pos <- "td"
  }
  if (isFALSE(cl.pos)) {
    cl.pos <- "n"
  }
  if (is.null(cl.pos) | isTRUE(cl.pos)) {
    if (type == "full") 
      cl.pos <- "r"
    if (type == "lower") 
      cl.pos <- "b"
    if (type == "upper") 
      cl.pos <- "r"
  }
  if (outline) 
    col.border <- "black"
  if (!outline) 
    col.border <- col.fill
  if (!add) {
    par(mar = mar, bg = "white")
    plot.new()
    xlabwidth <- ylabwidth <- 0
    for (i in 1:50) {
      xlim <- c(m1 - 0.5 - xlabwidth, m2 + 0.5 + mm * cl.ratio * 
                  (cl.pos == "r"))
      ylim <- c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b"), 
                n2 + 0.5 + ylabwidth)
      plot.window(xlim + c(-0.2, 0.2), ylim + c(-0.2, 0.2), 
                  asp = 1, xaxs = "i", yaxs = "i")
      x.tmp <- max(strwidth(newrownames, cex = tl.cex))
      y.tmp <- max(strwidth(newcolnames, cex = tl.cex))
      if (min(x.tmp - xlabwidth, y.tmp - ylabwidth) < 1e-04) 
        break
      xlabwidth <- x.tmp
      ylabwidth <- y.tmp
    }
    if (tl.pos == "n" | tl.pos == "d") 
      xlabwidth <- ylabwidth <- 0
    if (tl.pos == "td") 
      ylabwidth <- 0
    if (tl.pos == "ld") 
      xlabwidth <- 0
    laboffset <- strwidth("W", cex = tl.cex) * tl.offset
    xlim <- c(m1 - 0.5 - xlabwidth - laboffset, m2 + 0.5 + 
                mm * cl.ratio * (cl.pos == "r")) + c(-0.35, 0.15)
    ylim <- c(n1 - 0.5 - nn * cl.ratio * (cl.pos == "b"), 
              n2 + 0.5 + ylabwidth * abs(sin(tl.srt * pi/180)) + 
                laboffset) + c(-0.15, 0.35)
    if (.Platform$OS.type == "windows") {
      windows.options(width = 7, height = 7 * diff(ylim)/diff(xlim))
    }
    plot.window(xlim = xlim, ylim = ylim, asp = 1, xlab = "", 
                ylab = "", xaxs = "i", yaxs = "i")
  }
  laboffset <- strwidth("W", cex = tl.cex) * tl.offset
  symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
                                                         len.DAT), bg = bg, fg = bg)
  if (method == "circle" & plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, bg = col.fill, 
            circles = 0.9 * abs(DAT)^0.5/2, fg = col.border)
  }
  if (method == "ellipse" & plotCI == "n") {
    ell.dat <- function(rho, length = 99) {
      k <- seq(0, 2 * pi, length = length)
      x <- cos(k + acos(rho)/2)/2
      y <- cos(k - acos(rho)/2)/2
      return(cbind(rbind(x, y), c(NA, NA)))
    }
    ELL.dat <- lapply(DAT, ell.dat)
    ELL.dat2 <- 0.85 * matrix(unlist(ELL.dat), ncol = 2, 
                              byrow = TRUE)
    ELL.dat2 <- ELL.dat2 + Pos[rep(1:length(DAT), each = 100), 
                               ]
    polygon(ELL.dat2, border = col.border, col = col.fill)
  }
  if (method == "number" & plotCI == "n") {
    text(Pos[, 1], Pos[, 2], font = 2, col = col.fill, labels = round((DAT - 
                                                                         int) * ifelse(addCoefasPercent, 100, 1)/zoom, ifelse(addCoefasPercent, 
                                                                                                                              0, 2)))
  }
  if (method == "pie" & plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, circles = rep(0.5, 
                                                           len.DAT) * 0.85)
    pie.dat <- function(theta, length = 100) {
      k <- seq(pi/2, pi/2 - theta, length = 0.5 * length * 
                 abs(theta)/pi)
      x <- c(0, cos(k)/2, 0)
      y <- c(0, sin(k)/2, 0)
      return(cbind(rbind(x, y), c(NA, NA)))
    }
    PIE.dat <- lapply(DAT * 2 * pi, pie.dat)
    len.pie <- unlist(lapply(PIE.dat, length))/2
    PIE.dat2 <- 0.85 * matrix(unlist(PIE.dat), ncol = 2, 
                              byrow = TRUE)
    PIE.dat2 <- PIE.dat2 + Pos[rep(1:length(DAT), len.pie), 
                               ]
    polygon(PIE.dat2, border = "black", col = col.fill)
  }
  if (method == "shade" & plotCI == "n") {
    addshade <- match.arg(addshade)
    symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
                                                           len.DAT), bg = col.fill, fg = addgrid.col)
    shade.dat <- function(w) {
      x <- w[1]
      y <- w[2]
      rho <- w[3]
      x1 <- x - 0.5
      x2 <- x + 0.5
      y1 <- y - 0.5
      y2 <- y + 0.5
      dat <- NA
      if ((addshade == "positive" || addshade == "all") & 
          rho > 0) {
        dat <- cbind(c(x1, x1, x), c(y, y1, y1), c(x, 
                                                   x2, x2), c(y2, y2, y))
      }
      if ((addshade == "negative" || addshade == "all") & 
          rho < 0) {
        dat <- cbind(c(x1, x1, x), c(y, y2, y2), c(x, 
                                                   x2, x2), c(y1, y1, y))
      }
      return(t(dat))
    }
    pos_corr <- rbind(cbind(Pos, DAT))
    pos_corr2 <- split(pos_corr, 1:nrow(pos_corr))
    SHADE.dat <- matrix(na.omit(unlist(lapply(pos_corr2, 
                                              shade.dat))), byrow = TRUE, ncol = 4)
    segments(SHADE.dat[, 1], SHADE.dat[, 2], SHADE.dat[, 
                                                       3], SHADE.dat[, 4], col = shade.col, lwd = shade.lwd)
  }
  if (method == "square" & plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, squares = abs(DAT)^0.5, 
            bg = col.fill, fg = col.border)
  }
  if (method == "color" & plotCI == "n") {
    symbols(Pos, add = TRUE, inches = FALSE, squares = rep(1, 
                                                           len.DAT), bg = col.fill, fg = col.border)
  }
  symbols(Pos, add = TRUE, inches = FALSE, bg = NA, squares = rep(1, 
                                                                  len.DAT), fg = addgrid.col)
  if (plotCI != "n") {
    if (is.null(lowCI.mat) || is.null(uppCI.mat)) 
      stop("Need lowCI.mat and uppCI.mat!")
    if (!order == "original") {
      lowCI.mat <- lowCI.mat[ord, ord]
      uppCI.mat <- uppCI.mat[ord, ord]
    }
    pos.lowNew <- getPos.Dat(lowCI.mat)[[1]]
    lowNew <- getPos.Dat(lowCI.mat)[[2]]
    pos.uppNew <- getPos.Dat(uppCI.mat)[[1]]
    uppNew <- getPos.Dat(uppCI.mat)[[2]]
    if (!(method == "circle" || method == "square")) 
      stop("method shoud be circle or square if draw confidence interval!")
    k1 <- (abs(uppNew) > abs(lowNew))
    bigabs <- uppNew
    bigabs[which(!k1)] <- lowNew[!k1]
    smallabs <- lowNew
    smallabs[which(!k1)] <- uppNew[!k1]
    sig <- sign(uppNew * lowNew)
    if (plotCI == "circle") {
      symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
              inches = FALSE, circles = 0.95 * abs(bigabs)^0.5/2, 
              bg = ifelse(sig > 0, col.fill, col[ceiling((bigabs + 
                                                            1) * length(col)/2)]), fg = ifelse(sig > 0, 
                                                                                               col.fill, col[ceiling((bigabs + 1) * length(col)/2)]))
      symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
              inches = FALSE, circles = 0.95 * abs(smallabs)^0.5/2, 
              bg = ifelse(sig > 0, bg, col[ceiling((smallabs + 
                                                      1) * length(col)/2)]), fg = ifelse(sig > 0, 
                                                                                         col.fill, col[ceiling((smallabs + 1) * length(col)/2)]))
    }
    if (plotCI == "square") {
      symbols(pos.uppNew[, 1], pos.uppNew[, 2], add = TRUE, 
              inches = FALSE, squares = abs(bigabs)^0.5, bg = ifelse(sig > 
                                                                       0, col.fill, col[ceiling((bigabs + 1) * length(col)/2)]), 
              fg = ifelse(sig > 0, col.fill, col[ceiling((bigabs + 
                                                            1) * length(col)/2)]))
      symbols(pos.lowNew[, 1], pos.lowNew[, 2], add = TRUE, 
              inches = FALSE, squares = abs(smallabs)^0.5, 
              bg = ifelse(sig > 0, bg, col[ceiling((smallabs + 
                                                      1) * length(col)/2)]), fg = ifelse(sig > 0, 
                                                                                         col.fill, col[ceiling((smallabs + 1) * length(col)/2)]))
    }
    if (plotCI == "rect") {
      rect.width <- 0.25
      rect(pos.uppNew[, 1] - rect.width, pos.uppNew[, 2] + 
             smallabs/2, pos.uppNew[, 1] + rect.width, pos.uppNew[, 
                                                                  2] + bigabs/2, col = col.fill, border = col.fill)
      segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                                                        2] + DAT/2, pos.lowNew[, 1] + rect.width, pos.lowNew[, 
                                                                                                             2] + DAT/2, col = "black", lwd = 1)
      segments(pos.uppNew[, 1] - rect.width, pos.uppNew[, 
                                                        2] + uppNew/2, pos.uppNew[, 1] + rect.width, 
               pos.uppNew[, 2] + uppNew/2, col = "black", lwd = 1)
      segments(pos.lowNew[, 1] - rect.width, pos.lowNew[, 
                                                        2] + lowNew/2, pos.lowNew[, 1] + rect.width, 
               pos.lowNew[, 2] + lowNew/2, col = "black", lwd = 1)
      segments(pos.lowNew[, 1] - 0.5, pos.lowNew[, 2], 
               pos.lowNew[, 1] + 0.5, pos.lowNew[, 2], col = "grey70", 
               lty = 3)
    }
  }
  if (!is.null(p.mat) & !insig == "n") {
    if (!order == "original") 
      p.mat <- p.mat[ord, ord]
    pos.pNew <- getPos.Dat(p.mat)[[1]]
    pNew <- getPos.Dat(p.mat)[[2]]
    ind.p <- which(pNew > (sig.level))
    if (insig == "pch") {
      points(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
             pch = pch, col = pch.col, cex = pch.cex, lwd = 2)
    }
    if (insig == "p-value") {
      text(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
           round(pNew[ind.p], 2), col = pch.col)
    }
    if (insig == "blank") {
      symbols(pos.pNew[, 1][ind.p], pos.pNew[, 2][ind.p], 
              inches = FALSE, squares = rep(1, length(pos.pNew[, 
                                                               1][ind.p])), fg = addgrid.col, bg = bg, add = TRUE)
    }
  }
  if (cl.pos != "n") {
    colRange <- assign.color(cl.lim2)
    ind1 <- which(col == colRange[1])
    ind2 <- which(col == colRange[2])
    colbar <- col[ind1:ind2]
    if (is.null(cl.length)) 
      cl.length <- ifelse(length(colbar) > 20, 11, length(colbar) + 
                            1)
    labels <- seq(cl.lim[1], cl.lim[2], length = cl.length)
    at <- seq(0, 1, length = length(labels))
    if (cl.pos == "r") {
      vertical <- TRUE
      xlim <- c(m2 + 0.5 + mm * 0.02, m2 + 0.5 + mm * cl.ratio)
      ylim <- c(n1 - 0.5, n2 + 0.5)
    }
    if (cl.pos == "b") {
      vertical <- FALSE
      xlim <- c(m1 - 0.5, m2 + 0.5)
      ylim <- c(n1 - 0.5 - nn * cl.ratio, n1 - 0.5 - nn * 
                  0.02)
    }
    colorlegend(colbar = colbar, labels = round(labels, 2), 
                offset = cl.offset, ratio.colbar = 0.3, cex = cl.cex, 
                xlim = xlim, ylim = ylim, vertical = vertical, align = cl.align.text)
  }
  if (tl.pos != "n") {
    ylabwidth2 <- strwidth(newrownames, cex = tl.cex)
    xlabwidth2 <- strwidth(newcolnames, cex = tl.cex)
    pos.xlabel <- cbind(m1:m2, n2 + 0.5 + laboffset)
    pos.ylabel <- cbind(m1 - 0.5, n2:n1)
    if (tl.pos == "td") {
      if (type != "upper") 
        stop("type should be \"upper\" if tl.pos is \"dt\".")
      pos.ylabel <- cbind(m1:(m1 + nn) - 0.5, n2:n1)
    }
    if (tl.pos == "ld") {
      if (type != "lower") 
        stop("type should be \"lower\" if tl.pos is \"ld\".")
      pos.xlabel <- cbind(m1:m2, n2:(n2 - mm) + 0.5 + laboffset)
    }
    if (tl.pos == "d") {
      pos.ylabel <- cbind(m1:(m1 + nn) - 0.5, n2:n1)
      pos.ylabel <- pos.ylabel[1:min(n, m), ]
      symbols(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], add = TRUE, 
              bg = bg, fg = addgrid.col, inches = FALSE, squares = rep(1, 
                                                                       length(pos.ylabel[, 1])))
      text(pos.ylabel[, 1] + 0.5, pos.ylabel[, 2], newcolnames[1:min(n, 
                                                                     m)], col = tl.col.x, cex = tl.cex, ...)
    }
    else {
      text(pos.xlabel[, 1], pos.xlabel[, 2], newcolnames, 
           srt = tl.srt, adj = ifelse(tl.srt == 0, c(0.5, 
                                                     0), c(0, 0)), col = tl.col.x, cex = tl.cex, offset = tl.offset, 
           ...)
      text(pos.ylabel[, 1], pos.ylabel[, 2], newrownames, 
           col = tl.col.y, cex = tl.cex, pos = 2, offset = tl.offset, 
           ...)
    }
  }
  title(title, ...)
  if (!is.null(addCoef.col) & (!method == "number")) {
    text(Pos[, 1], Pos[, 2], col = addCoef.col, labels = round((DAT - 
                                                                  int) * ifelse(addCoefasPercent, 100, 1)/zoom, ifelse(addCoefasPercent, 
                                                                                                                       0, 2)))
  }
  if (type == "full" & plotCI == "n" & !is.null(addgrid.col)) 
    rect(m1 - 0.5, n1 - 0.5, m2 + 0.5, n2 + 0.5, border = addgrid.col)
  if (!is.null(addrect) & order == "hclust" & type == "full") {
    corrRect.hclust(corr, k = addrect, method = hclust.method, 
                    col = rect.col, lwd = rect.lwd)
  }
  invisible(corr)
}


plot_correlation_matrix <- function(gene.df, values, palette, filename, title) { 
  
  #gene.df <- byGRB_list[[32]][[1]]
  #z <- as.matrix(byGRB_list[[32]][[2]])
  z <- as.matrix(values)
  
  #gene.df <- gene.df[gene.df$ens %in% colnames(z),]
  genes_no_transcripts <- sapply(colnames(z), function(x) {
    unlist(str_split(x, ":"))[1]
  })
  gene.df <- gene.df[gene.df$ens %in% genes_no_transcripts,]
  
  # col order is (l. blue, black, green, yellow, red): 
  # one for GRB trgt, two for GRB byst, three for GWAS out, four for GWAS byst, five for GWAS trgt
  gene.df$col <- ifelse(gene.df$category=="GRB-trgt", palette[1], 
                        ifelse(gene.df$category=="GRB-byst", palette[2],
                               ifelse(gene.df$category=="GWAS-out", palette[3],
                                      ifelse(gene.df$category=="GWAS-byst", palette[4], palette[5]))))
  
  # reorder columns by genomic coordinate of genes
  gene.df$order <- 1:nrow(gene.df)
  m <- merge(data.frame(genes_no_transcripts, names(genes_no_transcripts)), 
             gene.df[,c("ens", "order", "col")], by=1)
  m <- m[order(m$order),]
  z <- z[, m$names.genes_no_transcripts.] 
  
  # special case is only one enh, adding an empty line
  if(is.null(dim(z))==T) {
    z <- rbind(z, rep(0, length(z)))
    rownames(z) <- c(rownames(values), "Ignore_me")
  }

  # transform gene IDs into hgnc if one exists, and leave transcript name 
  colnames(z) <- sapply(colnames(z), function(x) {
  
    gen_tr <- unlist(str_split(x, ":")) # splits into gene and transcript ID   
    y <- gene.df[gene.df[,1]==gen_tr[1], 2]
    if(str_length(y)>0) {
      return(paste(y, gen_tr[2], sep=":"))
    }
    else {
      return(x)
    }
  })
  
  png(filename, width=1000, height=500)
  
  #corrplot::corrplot(z, method="circle", title=names(byGRB_list)[32], is.corr=F, mar=c(0,0,1,0), tl.col=c(palette_autumn_leaves_on_wood, rep("black", 10)))
  corrplot_2col_labels(z, method="circle", title=title, is.corr=F, mar=c(0,0,1,0), tl.col.x=m$col, tl.col.y="black")
  legend("topright", legend = c("GRB-trgt", "GRB-byst", "GWAS-out", "GWAS-byst", "GWAS-trgt"), text.col = palette, bty="n", cex=1.2)
  dev.off()
}


# function that takes a vector of expression for an enhancer and a promoter (df with 2 columns)
# assigns a score to each position where enhancer has non-zero expression:
# one to highest expressed promoter, 2 to second highest, etc. 
# the lower the score the more right-handed the distribution of enh. is among promoters ordered by ascending expression (the better enh. expression correlates to high expression in promoters)
# comparing observed score against iter random positionings of enh. will give empirical P value
# returns empirical p and number of issues with enhancer expressed
empP_for_PE_pair <- function(expression_for_pair, iter=10000, plot=F, gene_name="Gene") {
  
  # one pair (E-P) 
  expression_for_pair <- expression_for_pair[order(expression_for_pair[,2]),]
  expression_for_pair$x <- nrow(expression_for_pair):1
  colnames(expression_for_pair) <- c("enh", "prom", "prom_order_asc")
  
  # which ones are E+ (have non-zero expression at enhancer),
  # make scoring function assigning 1, 2, 3, ... nrow(pair_vector) to first, second, ... highest expressed in promoter expression
  # this is the observed score for that enh-promoter pair (expected to be right-heavy)
  scoring_enhancers_from_highest <- function(vec) {
    vec$x <- nrow(vec):1
    return(sum(vec[,3][vec[,1]>0]))
  }
  
  # permute positions of E+ iter times to get distribution of scores
  how_many <- length(which(expression_for_pair$enh>0)) # how many tissues have non-zero expression of enhancer
  
  # makes iter random samplings of size how_many from enhancer, scores each
  distr <- integer(iter)
  for(i in 1:iter) {
    distr[i] <- sum(sample(1:nrow(expression_for_pair), size=how_many, replace=F))
    
  }
  
  # how many out of iter have lower score than observed -> empirical P
  p_emp <- length(which(distr<=scoring_enhancers_from_highest(expression_for_pair))) / iter
  
  if(plot==T) {
    
    plot_empP_for_pair(expression_for_pair, p_emp)
  }
  
  return(p_emp)
}


#first column is enhancer expression (data.frame), second is gene, 
#name of second column is used as geneID on the plot, first for enhancer
plot_empP_for_pair <- function(expression_for_pair, p_emp) {
  
  png(paste0(colnames(expression_for_pair)[2], "_", 
             colnames(expression_for_pair)[1], ".png"), 
      height=500, width=500)
  par(mar=c(1,2,2,1))
  layout(matrix(1:2, 2, 1), heights=c(4, 1))
  
  expression_for_pair <- expression_for_pair[order(expression_for_pair[,2]),]
  
  plot(1:nrow(expression_for_pair), expression_for_pair[,2], 
       type="l", xlab="FANTOM5 samples", ylab="Promoter expression", 
       main=colnames(expression_for_pair)[2], cex.main=1.5, xaxt='n')
  legend("topleft", legend=paste("P_emp:", p_emp, sep=" "), bty="n", cex=1.5)
  
  plot(0, xlim=c(1,nrow(expression_for_pair)), ylim=c(0,1), type='n', 
       main=colnames(expression_for_pair)[1], cex.main=1.5,
       xaxt='n', yaxt='n')
  abline(v=which(expression_for_pair[,1]>0), col="red")
  
  dev.off()
}



######################
# inputs GRanges with one entry per transcript
# each transcript is mapped to the highest scoring CAGE cluster (TPMs)
# will filter our all transcripts for the same gene with less then ratio (set to 1 percent) 
# of highest TPM cluster
filter_low_clusters <- function(gr, ratio=0.01) {
  
  gr$drop <- F
  
  genes_with_multiples <- unique(gr$ensembl_gene_id[which(duplicated(gr$ensembl_gene_id))])
  for (i in 1:length(genes_with_multiples)) {

    w <- which(gr$ensembl_gene_id==genes_with_multiples[i])
    m <- max(gr$tpm[w])
    gr$drop[w] <- ifelse(gr$tpm[w] >= m*ratio, F, T) # if tpm is <ratio*max.tpm, flags for dropping
  }
  gr <- gr[gr$drop==F]
  gr$drop <- NULL
  
  return(gr)
}


######################
# requires fantom5.enhancers.tpm.gr and promoter_cage_tpm_df preloaded
# make sure to use empP files from right folder (cageTSSs)
plots_empP_per_locus_from_scratch <- function(byGRBlist, i) {
  
  print(names(byGRBlist)[i])
  
  enh_vec <- rownames(byGRBlist[[i]][[2]])
  
  if(length(enh_vec)==0) {
    print ("No enhancer peaks found.")
  }
  else {
  
    # [[1]] is ordered by genomic coor, [[2]] isn't
    df1 <- data.frame(byGRBlist[[i]][[1]]$ens, 1:length(byGRBlist[[i]][[1]]$ens),
                     stringsAsFactors=F)
    df2 <- data.frame(sapply(colnames(byGRBlist[[i]][[2]]), function(x){ unlist(str_split(x, ":"))[1] }), 
                      colnames(byGRBlist[[i]][[2]]), stringsAsFactors=F)
    df <- merge(df1, df2, by=1, all=F)
  
    genes_ens <- df[order(df[,2]),1]
    genes_hgnc <- sapply(genes_ens, function(e) {
      byGRBlist[[i]][[1]]$symbol[byGRBlist[[i]][[1]]$ens==e]
    })

  
    # name is either gene (if there is only one transcript)
    # or gene:TSS1, gene:TSS2, ... if there are multiple transcripts for it
    df$name_fixed <- NA
    duplicates <- names(table(df[,1])[table(df[,1])>1])
    for (p in 1:nrow(df)) {
    #does gene have more than one TSS?
    if(df[p,1] %in% duplicates) { 
      if (str_length(genes_hgnc[genes_ens==df[p,1]][1])>0) { #more TSSs, hgnc exists
        df$name_fixed[p] <- paste(genes_hgnc[genes_ens==df[p,1]][1], 
                                  unlist(str_split(df[p,3], ":"))[2], sep=":")
      }
      else { # more TSSs, hgnc doesn't exist
        df$name_fixed[p] <- df[p,3]
      }
    }
    else {
      if (str_length(genes_hgnc[genes_ens==df[p,1]][1])>0) { # one TSS hgnc exists
        df$name_fixed[p] <- genes_hgnc[genes_ens==df[p,1]][1]
      }
      else { # one tss, no hgnc symbol
        df$name_fixed[p] <- df[p,1]
      }
    }
  }
  
    current_chr <- unlist(str_split(enh_vec[1], ":"))[1]
    
    # load empP file for the chr
    empP_file <- paste("~/Documents/SchizophreniaSNPs/empP_all_files/cageTSSs/",
                       "empP_CAGEenhPermissive_hg19Ens_6.4Mb_", 
                       current_chr, ".txt", sep="")
    empP_df <- read.delim(empP_file)
    
    fantom5.enhancers.tpm.gr <- fantom5.enhancers.tpm.gr[seqnames(fantom5.enhancers.tpm.gr)==current_chr]
    
    ####### this should be added to plotting function
    png(filename=paste0(names(byGRBlist)[i], "_corrplot_empP.png"), 
        width=140*length(genes_ens)+70, 
        height=56*length(enh_vec)+112)
    print(paste0("Plotting empP for ", names(byGRBlist)[i]))
    
    #how to scale
    par(mar=c(0.1,0.5,0.5,0.5))
    # +1 column for enhancer labels, +1 row for gene expression
    layout(matrix(data=seq(1, (length(genes_ens)+1)*(length(enh_vec)+1)), 
                  nrow=length(enh_vec)+1, ncol=length(genes_ens)+1, byrow=T), 
           heights=c(2, rep(1, length(enh_vec))), 
           widths=c(0.5, rep(1, length(genes_ens))))
    #col_vec <- brewer.pal(9, "YlOrRd")[c(4,6,8)]
    col_vec <- c("darkgoldenrod1", "darkorange1", "red")
    
    # 1) plot header line with gene expressions
    # first label
    plot(0, xlim=c(0,1), ylim=c(0,1), type='n', 
         xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
    text(0.5, 0.5, "Genes", cex=1.8, font=2)
    # now row with genes
    a <- sapply(df[order(df[,2]),3], function(g) {
      p <- promoter_cage_tpm_df[rownames(promoter_cage_tpm_df)==g,]
      plot(1:length(p), sort(p), type="l", axes=F, xlab="", ylab="")
      #title(main=genes_hgnc[genes_ens==g], cex.main=2)
      
      # label is ens if there is no hgnc
      legend("topleft", legend=df$name_fixed[df[,3]==g], bty='n', cex=1.5)

      box(which="plot", lty="solid")
      polygon(x=c(1:length(p), length(p), 0), y=c(sort(p), 0, 0), col="black")
    })
    
    # 2) plot a matrix of enhancer distributions:
    for (enh_num in 1:length(enh_vec)) {
      e <- enh_vec[enh_num]
      
      # plot enhancer label
      plot(0, xlim=c(0,1), ylim=c(0,1), type='n', 
           xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
      text(0.5, 0.5, paste0("Enh", enh_num), cex=1.8, font=2)
      
      #now plot a line for that enhancer against all genes:
      egr <- GRanges(e)
      s <- subsetByOverlaps(fantom5.enhancers.tpm.gr, egr)
      enh_expression <- apply(as.data.frame(mcols(s)), 2, sum) 
      
      c <- sapply(df[order(df[,2]),3], function(g) {
        
        p <- promoter_cage_tpm_df[rownames(promoter_cage_tpm_df)==g,]
        # enh expression vector needs to be reordered to match asc. gene's expression
        expression_for_pair <- as.data.frame(cbind(enh_expression, p))
        expression_for_pair <- expression_for_pair[order(expression_for_pair[,2]),] 
        
        empP <- empP_df$empP[empP_df$enh==e & empP_df$gene==g]
        stripes_col <- ifelse(empP>=0.01, col_vec[1], 
                              ifelse(empP>=0.001, col_vec[2], col_vec[3]))
        
        plot(0, xlim=c(1, nrow(expression_for_pair)), ylim=c(0,1), type='n', 
             xaxt='n', yaxt='n', xlab='', ylab='')
        if(g %in% rownames(promoter_cage_tpm_df)) {
          abline(v=which(expression_for_pair[,1]>0), col=stripes_col)
          text(floor(length(p)/2), 0.5, round(empP, 3), cex=3.0, font=2)
        }
      })
    }
    
    dev.off()
  }
}


######################
# simple wrapper to plot legend with boxes on top of empty plot
plot_box_legend <- function(string, col_vec, filename, title, horiz) {
  
  if(horiz==T) {
    h <- 250
    w <- 700
  }
  else {
    h <- 350
    w <- 400
  }
  
  png(filename=filename, width=w, height=h)
  plot(0, xlim=c(0,1), ylim=c(0,1), type='n', xaxt='n', yaxt='n', xlab='', ylab='')
  legend("center", legend=string, fill=col_vec, border="black", 
         bty='n', cex=2.0, text.font=2, title=title, horiz=horiz)
  dev.off()
}


######################
# requires objects mart, byGRBlist, loci.gr and mm10.GRBs.subset
# plots Gviz for each locus and puts it in a .png file
# labels of GRB and locus track can be turned off (others are irresponsice to the flag)
plot_locus_Gviz <- function(i, dir, labels=T, byGRB_list, 
                            loci.gr, snps.gr, grbs, mart) {
  
  # setting a nice pallete of Grey, Green and Purple that match
  Grey   <- "#7D949B"
  Green  <- "#9CA94E"
  Purple <- "#A74F7C" 
  
  if(nrow(byGRB_list[[i]][[2]])>0) {
    
    current_chr <- unlist(str_split(rownames(byGRB_list[[i]][[2]][1]), ":"))[1]
    s <- sapply(rownames(byGRB_list[[i]][[2]]), function(x){ as.numeric(unlist(str_split(x, "[:-]"))[2]) }) # starts of all enhancers
    e <- sapply(rownames(byGRB_list[[i]][[2]]), function(x){ as.numeric(unlist(str_split(x, "[:-]"))[3]) }) # ends of all enhancers
    range_of_enhancers <- GRanges(paste0(current_chr, ":", min(s), "-", max(e))) 
    
    current_loci <- subsetByOverlaps(loci.gr, range_of_enhancers)
    names(current_loci) <- current_loci$gwas_genes
    current_GRB <- subsetByOverlaps(grbs, current_loci)
    names(current_GRB) <- names(byGRB_list)[i]
    
    # plot ideogram and scale
    itrack <- IdeogramTrack(genome="hg19", chromosome=current_chr)
    Gtrack <- GenomeAxisTrack(name=as.character(current_chr))
    
    if(str_length(current_GRB$targets)==0) {
      grbtrack <- AnnotationTrack(current_GRB, name="GRB human:mouse", 
                                  fill=Grey, col=Grey)
    }else{
      grbtrack <- AnnotationTrack(current_GRB, name="GRB human:mouse", 
                                  id=current_GRB$targets, fontcolor.feature="black", 
                                  fill=Grey, col=Grey)
      
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
    snpTrSize <-  ifelse(length(current_SNPs)<4, 1, 
                         ifelse(length(current_SNPs)<10, 2,
                                ifelse(length(current_SNPs)<20, 3,
                                       ifelse(length(current_SNPs)<30, 4, 
                                              ifelse(length(current_SNPs)<40, 5, 6)))))
    
    
    if(str_length(names(current_loci))==0) {
      ltrack <- AnnotationTrack(current_loci, name="GWAS locus",
                                fill=Purple, col=Purple)
      #fill=Grey, col=Grey)
    }
    else {
      ltrack <- AnnotationTrack(current_loci, name="GWAS locus", 
                                fill=Purple, col=Purple,
                                #fill=Grey, col=Grey,
                                id=names(current_loci), fontcolor.feature="black")
      #group=names(current_loci), groupAnnotation="group") # this throws a single genome warning - ignore!
      #displayPars(ltrack) <- list(showId=T)
      if(labels) {
        displayPars(ltrack) <- list(featureAnnotation="id")
      }
    }
    
    strack <- AnnotationTrack(current_SNPs, name="GWAS SNPs", 
                              group=current_SNPs$name, groupAnnotation="group",
                              fill=Purple, col=Purple) 
    if(labels) {
      displayPars(strack) <- list(showId=T)
    }
    
    # plot enhancer track
    splt <- str_split(rownames(byGRB_list[[i]][[2]]), "[:-]")
    e <- GRanges(seqnames=sapply(splt, function(x){ x[1] }),
                 ranges=IRanges(start=as.integer(sapply(splt, function(x) { x[2] })), 
                                end=as.integer(sapply(splt, function(x){ x[3] }))
                 ))
    etrack <- AnnotationTrack(e, name="CAGE enhancers", 
                              stacking="squish", # default, used to be 'dense', but even in 'full' it's not expanded
                              group=paste0("Enh", 1:length(e)), groupAnnotation="group",
                              #fill=Grey, col=Grey) # pink is brewer.pal(8, "PiYG")[1]
                              fill=Purple, col=Purple) # pink is brewer.pal(8, "PiYG")[1]
    if(labels) {
      displayPars(etrack) <- list(showId=T) # turn on group annotation
    }
    
    # resizing track with for enhancer track
    enhTrSize <-  ifelse(length(e)<4, 1, 
                         ifelse(length(e)<7, 2,
                                ifelse(length(e)<10, 3,
                                       ifelse(length(e)<16, 4, 
                                              ifelse(length(e)<25, 6, 8)))))
  
    # plot gene track 
    splt <- str_split(colnames(byGRB_list[[i]][[2]]), ":")
    tr_in_locus <- sapply(splt, function(x) { return(x[2]) }) # transcripts to keep for plotting
    
    gtrack <- BiomartGeneRegionTrack(start=start(current_GRB), end=end(current_GRB), chromosome=current_chr, 
                                     biomart=mart, name="Genes")
    gtrack <- gtrack[transcript(gtrack) %in% tr_in_locus]
    
    # edition: added unlist here
    feature(gtrack) <- unlist(sapply(gene(gtrack), function(g) { 
      byGRB_list[[i]][[1]]$category[byGRB_list[[i]][[1]]$ens==g]
    }))
    feature(gtrack) <- str_replace_all(feature(gtrack), "-", "")
    
    # resizing track with for genes
    gTrSize <-  ifelse(length(gtrack)<3, 1, 
                         ifelse(length(gtrack)<5, 2,
                                ifelse(length(gtrack)<8, 3,
                                       ifelse(length(gtrack)<16, 5, 
                                              ifelse(length(gtrack)<25, 7, 
                                                     ifelse(length(gtrack)<37, 9, 11))))))
    
    # using palette_primary_desktop <- c("#4E7CD9", "#153259", "#285911", "#D9A60D", "#D93425") here...
    col_palette <- palette_primary_desktop
    displayPars(gtrack) <- list(showId=T, col=NULL,
                                GWASout=col_palette[3], GWAStrgt=col_palette[5], GRBtrgt=col_palette[1], 
                                GRBbyst=col_palette[2], GWASbyst=col_palette[4]) # turn on gene symbols, no dashes allowed
    
    # plot CNE density
    CNEdensfile <- "~/Documents/grbs/HCNE_density_hg19_mm10_96pc_50col.wig.gz"
    CNEdtrack <- DataTrack(range=CNEdensfile, name="hg19_mm10\n96/50",
                           genome="hg19", chromosome=current_chr, 
                           type="polygon", col.mountain="azure3",
                           fill.mountain=c("azure3", "azure3"))
    
    # plot directionality index
    dirfile <- "~/Documents/snpsInGrbs/precalcData/tad_tracks/directionality_homer_neural.bigWig"
    dirtrack <- DataTrack(range=dirfile, name="Homer_DI_neural",
                          genome="hg19", chromosome=current_chr, 
                          type="mountain", baseline=0, lty.mountain=0,
                          fill.mountain=adjustcolor(c("darkgreen", "darkred"), alpha.f = 0.7),
                          col.mountain="black") # baseline line colour
    
    dirTrSize <- 1.5 # size of the dir. index track
    # high values look noisy if not smoothed, low ones throw an error if smoothed 
    # so only smooth of directionality index is above this arbitrary threshold
    DI_neural <- import.bw(dirfile)
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
    plotTracks(list(itrack, Gtrack, 
                    grbtrack, CNEdtrack, dirtrack,
                    ltrack, strack, etrack, gtrack), 
               from=xmin, to=xmax, 
               sizes=c(1, 1, 
                       1, 1.5, dirTrSize, 
                       1, snpTrSize, enhTrSize, gTrSize),
               background.title="azure3", fontcolor.title="black")
    dev.off()
  }
  else {
    print(paste("No Gviz plot for", names(byGRB_list)[i], "no enhancers found.", sep=" "))
  }
}
# how to collapse enhancers?
# add GRBs over 7 species


######################
plot_beanplot_individual <- function(expression_for_pair, neuronal, col2plot, sign="mwu", n) {
  
  par(mar=c(2,1,2,1))
  
  # significance can be determined by Mann-Whitney-U test, 
  # or permutation of the median based (needs raw, pre-log TPMs) 
  if(sign=="mwu") {
    # two-sided unpaired two sample wilcox test is MW-U test
    p <- wilcox.test(expression_for_pair$p[expression_for_pair$enh_expression==0], 
                     expression_for_pair$p[expression_for_pair$enh_expression> 0])$p.value
  }
  else if(sign=="permutation") {
    #determining significance string to plot (permutation-based deviation of the mean)
    p <- empP_diff_in_median_sampling(expression_for_pair, iter=10000)
  }
  else {
    p <- NA
  }
  # multiple testing correction, Bonferonni
  p <- n*p
  
  # ** in p <=0.01, *** if p<=0.001, nothing otherwise
  s <- ifelse(p<=0.001, "***", ifelse(p<=0.01, "**", ifelse(p<=0.05, "*", "")))
  
  expression_for_pair$enh_expression <- log_transform(expression_for_pair$enh_expression, 10)
  expression_for_pair$p <- log_transform(expression_for_pair$p, 10)
  
  beanplot(expression_for_pair$p[expression_for_pair$enh_expression==0], 
           expression_for_pair$p[expression_for_pair$enh_expression> 0],
           # default bw reports 'Error in bw.SJ(x, method = "dpi") : sample is too sparse to find TD'
           show.names=F, bw="nrd0", 
           #wd=bp$wd*(vec/max(vec)),
           log="", side="both",
           overallline = "median",
           beanlines="median",
           col= list(brewer.pal(8, "Dark2")[8], c(col2plot, "white")),
           what=c(0,1,1,0), border=NA, ylim=c(0,3)) 
  
  legend(x=0.75, y=3.75, legend=s, bty="n", cex=6, text.col="black")
}


######################
# subsample size is equal to fraction of tissues with non-zero expression of enhancer (how_many)
# creates distr. of medians, based on iter subsamples of the same size (how_many)
# returns single-sided empirical p: 
# a fraction of iterations with median greater or equal median observed for enhancer
empP_diff_in_median_sampling <- function(expression_for_pair, iter=10000) {
  
  # how many tissues have non-zero expression of enhancer
  how_many <- length(which(expression_for_pair$enh_expression>0)) 
  
  # makes iter random samplings of size how_many from enhancer, scores each
  distr <- integer(iter)
  for(i in 1:iter) {
    distr[i] <- median(sample(expression_for_pair$p, size=how_many, replace=F))
  }
  
  # how many out of iter have lower score than observed -> empirical P
  p_emp <- length(which(distr>=median(expression_for_pair$p[expression_for_pair$enh_expression>0]))) / iter
  return(p_emp)
}


######################
# requires promoter_cage_tpm_df and fantom5.enhancers.tpm.gr preloaded
# 
plots_prom_expression_beanplot_per_locus <- function(byGRBlist, i, neuronal=F, sign="permutation",
                                                     dir, file_suffix="_gene_expression_beanplot_mw.png") {
  
  print(names(byGRBlist)[i])
  enh_vec <- rownames(byGRBlist[[i]][[2]])
  
  if(length(enh_vec)==0) {
    print ("No enhancer peaks found.")
  }
  else {
    
    # [[1]] is ordered by genomic coor, [[2]] isn't
    df1 <- data.frame(byGRBlist[[i]][[1]]$ens, 1:length(byGRBlist[[i]][[1]]$ens),
                      stringsAsFactors=F)
    df2 <- data.frame(sapply(colnames(byGRBlist[[i]][[2]]), function(x){ unlist(str_split(x, ":"))[1] }), 
                      colnames(byGRBlist[[i]][[2]]), stringsAsFactors=F)
    df <- merge(df1, df2, by=1, all=F)
    
    genes_ens <- df[order(df[,2]),1]
    genes_hgnc <- sapply(genes_ens, function(e) {
      byGRBlist[[i]][[1]]$symbol[byGRBlist[[i]][[1]]$ens==e]
    })
    # total number of tests per locus, to be used for Bonferonni
    # correction in plot_beanplot_individual
    n <- length(genes_ens)*length(enh_vec) 
    
    # name is either gene (if there is only one transcript)
    # or gene:TSS1, gene:TSS2, ... if there are multiple transcripts for it
    df$name_fixed <- NA
    duplicates <- names(table(df[,1])[table(df[,1])>1])
    for (p in 1:nrow(df)) {
      #does gene have more than one TSS?
      if(df[p,1] %in% duplicates) { 
        if (str_length(genes_hgnc[genes_ens==df[p,1]][1])>0) { #more TSSs, hgnc exists
          df$name_fixed[p] <- paste(genes_hgnc[genes_ens==df[p,1]][1], 
                                    unlist(str_split(df[p,3], ":"))[2], sep=":")
        }
        else { # more TSSs, hgnc doesn't exist
          df$name_fixed[p] <- df[p,3]
        }
      }
      else {
        if (str_length(genes_hgnc[genes_ens==df[p,1]][1])>0) { # one TSS hgnc exists
          df$name_fixed[p] <- genes_hgnc[genes_ens==df[p,1]][1]
        }
        else { # one tss, no hgnc symbol
          df$name_fixed[p] <- df[p,1]
        }
      }
    }
    
    # there are actually some filenames longer than 255 chars
    # f <- substr(names(byGRBlist)[i], 1, (255-str_length(file_suffix)))
    # f <- paste0(dir, f, file_suffix)
    png(filename=paste0(dir, names(byGRBlist)[i], file_suffix), 
        width=200*length(genes_ens)+140, 
        height=200*length(enh_vec)+100)
    # consider
    #    png(filename=paste0(dir, names(byGRBlist)[i], file_suffix), 
   # width=300*length(genes_ens)+180, 
  #  height=200*length(enh_vec)+100)
    
    print(paste0("Plotting beanplot for ", names(byGRBlist)[i]))
    
    #how to scale matrix of beanplots
    par(mar=c(0.1,0.5,0.5,0.5))
    # +1 column for enhancer labels, +1 row for gene label
    layout(matrix(data=seq(1, (length(genes_ens)+1)*(length(enh_vec)+1)), 
                  nrow=length(enh_vec)+1, ncol=length(genes_ens)+1, byrow=T), 
           heights=c(0.5, rep(1, length(enh_vec))), 
           widths=c(0.7, rep(1, length(genes_ens)))) # change accordingly to layout up there
    
    # 1) just lists gene IDs
    plot(0, xlim=c(0,1), ylim=c(0,1), type='n', 
         xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
    text(0.5, 0.5, "Genes", cex=1.8, font=2)
    
    a <- sapply(df[order(df[,2]),3], function(g) {
      plot(0, xlim=c(0,1), ylim=c(0,1), type='n', 
           xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
      text(0.5, 0.5, df[df[,3]==g,4], cex=1.8, font=2)
    })
    
    # 2) plot a matrix of beanplots:
    for (enh_num in 1:length(enh_vec)) {
      
      e <- enh_vec[enh_num]
      s <- subsetByOverlaps(fantom5.enhancers.tpm.gr, GRanges(e))
      enh_expression <- apply(as.data.frame(mcols(s)), 2, sum) 
      
      # prepare counts of grey/coloured tissues
      if(neuronal==T) {
        vec <- c(length(enh_expression[names(enh_expression) %in% neuronal_samples$ids 
                                       & enh_expression==0]),
                 length(enh_expression[names(enh_expression) %in% neuronal_samples$ids 
                                       & enh_expression>0]))
        
        enh_expression <- enh_expression[names(enh_expression) %in% neuronal_samples$ids]
      }
      else {
        vec <- c(length(which(enh_expression==0)),
                 length(which(enh_expression>0)))
      }
      
      col2plot <- ifelse(neuronal==F,
                         brewer.pal(8, "PiYG")[1], # pink
                         brewer.pal(8, "PiYG")[8]) # green
      
      # plot enhancer label that has number of grey and pink points
      plot(0, xlim=c(0,1), ylim=c(0,1), type='n', 
           xaxt='n', yaxt='n', xlab='', ylab='', bty='n')
      text(0.5, 0.8, paste0("Enh", enh_num), cex=3.5, font=2)
      text(0.5, 0.5, paste("E0:", vec[1], sep=" "), cex=3.5, font=1.5, 
           col=brewer.pal(8, "Dark2")[8]) # dark gray
      text(0.5, 0.3, paste("E+:", vec[2], sep=" "), cex=3.5, font=1.5, 
           col=col2plot)
      
       # now plot beanplot
        a <- sapply(df[order(df[,2]),3], function(g) {
        
        p <- promoter_cage_tpm_df[rownames(promoter_cage_tpm_df)==g,]
        # will log-transform data and plot, and add significance based on sign method
        expression_for_pair <- data.frame(enh_expression=enh_expression, 
                                          p=p)
        plot_beanplot_individual(expression_for_pair, neuronal, col2plot, sign=sign, n=n)
      })
    }
    dev.off()
  }
}


######################
# adds one so the min. value is zero again
log_transform <- function(x, base) {
  x <- x+1
  return(log(x, base))
}


######################
# reverse of log_transform above
rev_log_transform <- function(logx, base) {
  x <- base^logx
  return(x-1)
}

######################
# returns a list where the first element is 
#every combination of enhancer in LD and transcript promoter in the GRB
# with E-TSS distance
# the second element is every enhancer-LD block combination and the trait from that LD
per_GRB_enhancers_classification <- function(enhancers.gr, current_GRB, genes_to_include, loci.gr, TSSs.perTranscript.gr) {
  
  l <- subsetByOverlaps(loci.gr, current_GRB)
  e <- granges(subsetByOverlaps(enhancers.gr, l))
  tr <- TSSs.perTranscript.gr[TSSs.perTranscript.gr$ensembl_gene_id %in% genes_to_include]
  e_l_t <- data.frame()
  e_p <- data.frame()
  
  if(length(e)>0) {
    #enhancer-LD-trait df
    e_l_t <- as.data.frame(findOverlaps(e, l))
    e_l_t$e <- as.character(e[e_l_t$queryHits])
    e_l_t$l <- as.character(l[e_l_t$subjectHits])
    e_l_t$t <- l$trait[e_l_t$subjectHits]
    e_l_t$queryHits <- NULL
    e_l_t$subjectHits <- NULL
    
    if(length(tr)>0){
      # enhancer-promoter distance df
      e_p <- expand.grid(as.character(e), names(tr))
      colnames(e_p) <- c("e", "tr")
      e_p$dist <- apply(e_p, 1, function(x){
        d <- distance(GRanges(x[1]), tr[grep(x[2], names(tr))])
      })
    }
  }
  return(list(e_p, e_l_t))
}



######################
# in FANTOM5 samples, is merging fetal and newborn legit? test and then merge....
# parietal lobe two [63:64] are 0.9 sp.corr
# spinal cord two [79:80] are 0.91 sp.corr
# renaming them so they appear as replicates
newborn_to_fetal <- function(nms){
  nms <- str_replace_all(nms, "newborn", "fetal")
}


######################
# used to pool FANTOM replicates/samples, 
# use with strng <- c("donor", "rep", "pool") 
# merging samples by splitting nms by strng (this is order-dependent), 
# and then merging all with same first element
# order in strng is important!
  merge_list_elements_before_word <- function(nms, strng) {
    biolRepList <- list()
    
    nms <- newborn_to_fetal(nms)
    
    # sort out capitalization and remove treated CD14 cells
    other <- sort_out_fantom5_filenames(nms)
    
    for(s in strng){
      t <- grep(s, other, value=T)
      other <- grep(s, other, invert=T, value=T)
      
      for(tis in unique(sapply(t, function(x){ str_split(x, s, simplify=T)[1] }))){
        ts <- grep(tis, t, value=T)
        biolRepList[[tis]] <- ts
      }
    } 
    # here other has all nms not caught by words in nms, presumably with no replicates
    
    l <- list()
    l <- c(l, other)
    names(l) <- other
    biolRepList <- c(biolRepList, l)
    return(biolRepList)
  }


######################
# used when merging FANTOM5 samples...
sort_out_fantom5_filenames <- function(nms){
  nms <- str_replace_all(nms, " - ", " ")
  nms <- str_replace_all(nms, "%2c ", " ")
  nms <- str_trim(nms)
  
  # there have been some inconsistent naming (upper/lower case)
  nms <- tolower(nms)
  nms <- str_replace_all(nms, "cnhs", "CNhs")
  
  # to remove bunch of similar CD14+ samples upon different treatments
  nms <- grep("treated", nms, 
              value=T, ignore.case=T, invert=T)
  
  return(nms)
}

######################
# from the full FANTOM5 name, extracts just the "CNhs" bit
get_FANTOM5_ids <- function(nms){
  ids <- sapply(nms, function(nm){
    s <- str_split(nm, "CNhs", simplify=T)
    s2 <- str_split(s[2], "\\.", simplify=T)[1]
    paste0("CNhs", s2)
  })
  return(ids)
}


######################
# all values here are already log-transformed, so it's the geom. mean of TPMs that I'm getting here
# need to revert to get the promoter_cage_tpm_df
merge_replicate_log_tpms <- function(l, lp) {
  
  # 0 cols df
  df <- data.frame(matrix(vector(mode='numeric', length=0), nrow=nrow(lp), ncol=0))
  rownames(df) <- rownames(lp)
  
  for(i in 1:length(l)){
    print(names(l)[i])
    ids <- get_FANTOM5_ids(l[[i]])
    
    if(length(ids)>1 & any(ids!="CNhsNA")){
      
      # calculates mean log expression and adds a column with it
      d <- rowMeans(lp[,ids])
      df <- cbind(df, d)
      
      # adds colname with : separated all ids
      colnames(df)[ncol(df)] <- paste(ids, collapse=":")
    }
    else{ # adds the sample with no replicates
      
      if(ids!="CNhsNA"){
        df <- cbind(df, lp[,ids])
        colnames(df)[ncol(df)] <- ids
      }
    }
  }
  
  return(df)
}


