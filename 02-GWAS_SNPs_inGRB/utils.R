# library(rtracklayer)
# library(stringr)
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(corrplot)
# library(plyr)
# library(Gviz)
# library(RColorBrewer)
# library(beanplot)
# options(scipen=999)

assemblies <- c("canFam3", "mm10", "monDom5", "galGal4", "lepOcu1", "danRer7", "tetNig2")

palette_primary_desktop <- c("#4E7CD9", "#153259", "#285911", "#D9A60D", "#D93425")

# Making granges from dataframe
grangesFromDataframe <- function(data, colName, keep_metadata=TRUE){
  # Description: make granges from data frame
  # which has a column of this format:
  # chr4:11394774-11431389
    # data = gRanges object
    # colName = column which contains the range
    # keep_metadata = store metadata in the dataframe
  
  data_r <-
    data %>%
    mutate(start=as.character(str_split(!!rlang::enquo(colName), pattern = ":|-", simplify = TRUE)[,2]),
           end=as.character(str_split(!!rlang::enquo(colName), pattern = ":|-", simplify = TRUE)[,3]),
           CHR= as.character(str_split(!!rlang::enquo(colName), pattern = ":|-", simplify = TRUE)[,1]))
  
  
  # Granges object is needed to detect SNP loci that overlap 2 GRBs
  data_r <-
    makeGRangesFromDataFrame(data_r,
                             seqnames.field ="chr",
                             start.field = "start",
                             end.field = "end",
                             keep.extra.columns = keep_metadata)
  
  
  return(data_r)
}


# current_GRB is a GRange of current GRB to be processed
# current_loci needs to have a GRanges with ranges of the loci overlapping that GRB
# and (colon delimited) string of GWAS associated in $gwas_genes field, and (colon delimited) string of target genes in $GRBname field (all target s predicted in all 7 pairwise comparisons)
# outputs a df with 3 columns, one row for each gene in that GRB, 
# plus any GWAS-associated genes with loci overlapping that GRB that fall outside of the GRB
# first column is EnsemblID, second id hgnc, third can be one of 5 classes
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
