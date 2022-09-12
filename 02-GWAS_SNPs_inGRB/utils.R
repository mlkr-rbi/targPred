# library(rtracklayer)
# library(stringr)
# library(BSgenome.Hsapiens.UCSC.hg19)
# library(corrplot)
# library(plyr)
# library(Gviz)
library(RColorBrewer)
library(beanplot)
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


load_per_chr_empP_file <- function(current_chr) {
  
  f <- paste0("data/empP_files/empP_CAGEenhPermissive_hg19Ens_6.4Mb_", as.character(current_chr), ".txt")
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


### ### ###

# requires promoter_cage_tpm_df and fantom5.enhancers.tpm.gr preloaded
# 
plots_prom_expression_beanplot_per_locus <- function(byGRBlist, i, neuronal=F, sign="permutation",
                                                     dir, file_suffix="_gene_expression_beanplot_mw.png") {
  
  print(names(byGRBlist)[i])
  enh_vec <- rownames(byGRBlist[[i]][[3]])
  
  if(length(enh_vec)==0) {
    print ("No enhancer peaks found.")
  }
  else {
    
    # [[1]] is ordered by genomic coor, [[2]] isn't
    df1 <- data.frame(byGRBlist[[i]][[1]]$ens, 1:length(byGRBlist[[i]][[1]]$ens),
                      stringsAsFactors=F)
    df2 <- data.frame(sapply(colnames(byGRBlist[[i]][[3]]), function(x){ unlist(str_split(x, ":"))[1] }), 
                      colnames(byGRBlist[[i]][[3]]), stringsAsFactors=F)
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

log_transform <- function(x, base) {
  x <- x+1
  return(log(x, base))
}
