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


fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}

give.n <- function(x){
  return(c(y = -0.5, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

ProportionOfGenesShared <- function(set1, set2){
  union_set <- union(set1,  set2)
  inters_set <- intersect(set1, set2)
  if(any(is.element(inters_set, union_set))){
    # return(table(is.element(inters_set, union_set))[[TRUE]]/length(union_set))}
    return(length(inters_set)/length(union_set))}
  else{return(0)}
  
}

GeneSetProportionMatrix <- function(x, exclude_set=0){
  # First col contains phenotype
  # Second col contains comma separated genes
  n <- nrow(x)
  len_vec <- c()
  
  gene_vec <- list()
  for(i in 1:n){
    gene_vec <- append(gene_vec, str_split(x[i,2], pattern = ", "))
    len_vec <- append(len_vec, length(gene_vec[[length(gene_vec)]]))
  }
  prop_matx <- matrix(nrow = n, ncol = n)
  for(i in c(1:n)){
    for(j in c(1:n)){
      prop_matx[i,j] <- ProportionOfGenesShared(gene_vec[[i]], gene_vec[[j]])
    }
  }
  # print(gene_vec)
  # print(len_vec)
  rownames(prop_matx) <- data.frame(x)[,1]
  
  prop_matx <- prop_matx[(len_vec > exclude_set),(len_vec > exclude_set) ]
  return(prop_matx)
}


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
    mutate(start=as.character(str_split(!!rlang::enquo(colName), pattern = ":|-|_", simplify = TRUE)[,2]),
           end=as.character(str_split(!!rlang::enquo(colName), pattern = ":|-|_", simplify = TRUE)[,3]),
           CHR= as.character(str_split(!!rlang::enquo(colName), pattern = ":|-|_", simplify = TRUE)[,1]))
  
  
  # Granges object is needed to detect SNP loci that overlap 2 GRBs
  data_r <-
    makeGRangesFromDataFrame(data_r,
                             seqnames.field ="chr",
                             start.field = "start",
                             end.field = "end",
                             keep.extra.columns = keep_metadata)
  
  
  return(data_r)
}

# Make BED file from granges with a name column
makeBed <- function(gr, name = FALSE, filename, 
                    includetrackname = FALSE, 
                    track_name ="My Track", 
                    track_description = "My description",
                    # since bed is 0 based and granges is 1 based
                    modify_start = 0){
  
  
  if(length(name) == 1 && name == FALSE){
    name <- c(rep(".", length(gr)))
  }
  else if(any(name=="")){
    message("Empty strings are not allowed, replacing them with NA")
    name <- base::replace(name, name == "", NA)
  }
  
  df <- data.frame(seqnames=seqnames(gr),
                   # Do not need to lessen it by one since when importing bed file it was increased
                   # starts=start(grb)-1,
                   starts=start(gr) + modify_start,
                   ends=end(gr),
                   names=name,
                   scores=c(rep(0, length(gr))),
                   strands=c(rep(".", length(gr))))
  
  write.table(df, quote=F, sep="\t", row.names=F, col.names=F, 
              file = filename)
  
  if(includetrackname){
    track_full = paste("track name=\"", track_name, "\" description=\"", track_description, "\" itemRgb=\"On\"")
    con <- file(filename, open="wt")
    writeLines(track_full, con)
    write.table(df, quote=F, sep="\t", row.names=F, col.names=F, file = con)
    close(con)
      }
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
  target <- data.frame(symbol=unlist(str_split(current_loci$GRB_GENES, ",")), 
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


load_per_chr_empP_file <- function(current_chr, emp_dir = "", fn="empP_CAGEenhPermissive_hg19Ens_6.4Mb_" ) {
  
  f <- file.path(emp_dir, paste0(as.character(current_chr), fn,  ".txt"))
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

writeGrangesToText <- function(gr, dir = "", is_empP_file = FALSE){
  
  if(class(gr) != "GRanges"){
    errorCondition("The object is not of class GRanges")
  }
  
  if(is_empP_file){
    coor <- "enh"
  }
  else{coor <- "coordinates"}
  mcol_names <- colnames(mcols(gr))
  
  gr.df <- as.data.frame(gr) 
  gr.df[coor] <- paste0(gr.df$seqnames, ":", gr.df$start, "-", gr.df$end)
  gr.df$chr <- gr.df$seqnames
  # Reordering columns 
  gr.df <- gr.df[, c("chr", coor, mcol_names)]
  write.table(gr.df, file = dir, quote = FALSE, sep = "\t", row.names = FALSE)
  
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
plots_prom_expression_beanplot_per_locus <- 
  function(byGRBlist,
           i, 
           neuronal=F, 
           sign="permutation",
           dir,
           file_suffix="_gene_expression_beanplot_mw.png") {
    
  # description: add description and paramater definitions
    
  enh_vec <- rownames(byGRBlist[[i]][[3]])
  
  if(length(enh_vec)==0) {
    print ("No enhancer peaks found.")
  }
  # Maybe remove this else and add a loop break in if
  else {
    
    # [[1]] is ordered by genomic coor, [[2]] isn't
    # Preparing df
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
    
    # Preparing plot
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
      print(paste("E0:", vec[1], sep=" "))
      print(paste("E+:", vec[2], sep=" ", enh_num))
      # now plot beanplot
      a <- sapply(df[order(df[,2]),3], function(g) {
        print(paste("enh: ", e, " prom ", g))
        p <- promoter_cage_tpm_df[rownames(promoter_cage_tpm_df)==g,]
        # will log-transform data and plot, and add significance based on sign method
        expression_for_pair <- data.frame(enh_expression=enh_expression, 
                                          prom=p)
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
    p <- wilcox.test(expression_for_pair$prom[expression_for_pair$enh_expression==0], 
                     expression_for_pair$prom[expression_for_pair$enh_expression> 0])$p.value
  }
  else if(sign=="permutation") {
    #determining significance string to plot (permutation-based deviation of the mean)
    p <- empP_diff_in_median_sampling(expression_for_pair, iter=10000)
    print(p)
  }
  else {
    p <- NA
  }
  # multiple testing correction, Bonferonni
  p <- n*p
  
  # ** in p <=0.01, *** if p<=0.001, nothing otherwise
  s <- ifelse(p<=0.001, "***", ifelse(p<=0.01, "**", ifelse(p<=0.05, "*", "")))
  #print(paste("P after correction", p))
  print(paste("Dimension of expression pair", dim(expression_for_pair)))
  
  expression_for_pair$enh_expression <- log_transform(expression_for_pair$enh_expression, 10)
  expression_for_pair$prom <- log_transform(expression_for_pair$prom, 10)
  
  beanplot(expression_for_pair$prom[expression_for_pair$enh_expression==0], 
           expression_for_pair$prom[expression_for_pair$enh_expression> 0],
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


# requires promoter_cage_tpm_df and fantom5.enhancers.tpm.gr preloaded
plotEnhPromBeanPerTrait <- 
  function(byGRBlist,
           i, 
           neuronal=F, 
           sign="permutation",
           loci.gr,
           dir,
           file_suffix="_gene_expression_beanplot.png",
           enhancer_tpm.gr,
           promoter_tpm.df) {
    
    # description: add description and parameter definitions
    
    # Vector of enhancers in a GRB
    enh_vec <- rownames(byGRBlist[[i]][[3]])
    #return(print(enh_vec))
    # Generating granges of enhancers in order to find overlaps with trait loci
    if(length(enh_vec) != 0){
    #print("No enhancer peaks found.")

    enh_vec.gr <- grangesFromDataframe(as.data.frame(enh_vec), colName = enh_vec)
    # Vector of enhancers that overlap trait associated loci
    enh_vec.gr <- subsetByOverlaps(enh_vec.gr, loci.gr)
    
    enh_vec <- mcols(enh_vec.gr)$enh_vec
    print(enh_vec)
    }
    

    if(length(enh_vec)==0) {
      print ("No enhancer peaks found.")
    }
    # Maybe remove this else and add a loop break in if
    else {
      
      # [[1]] is ordered by genomic coor, [[2]] isn't
      # Preparing df
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
      # correction in beanplotIndividual
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
      
      # Preparing plot
      # there are actually some filenames longer than 255 chars
      # f <- substr(names(byGRBlist)[i], 1, (255-str_length(file_suffix)))
      # f <- paste0(dir, f, file_suffix)
      png(filename=paste0(dir,"/",file_suffix), 
          width=200*length(genes_ens)+140, 
          height=200*length(enh_vec)+100)
      # consider
      #    png(filename=paste0(dir, names(byGRBlist)[i], file_suffix), 
      # width=300*length(genes_ens)+180, 
      #  height=200*length(enh_vec)+100)
      
      print(paste0("Plotting beanplot for ", i))
      
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
      #start_time <- Sys.time()
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
          
          prom_expression <- promoter_tpm.df[rownames(promoter_tpm.df)==g,]
          # will log-transform data and plot, and add significance based on sign method
          expression_for_pair <- data.frame(enh_expression=enh_expression, 
                                            prom=prom_expression)

          p_val <- byGRBlist[[i]][[3]][e,g]
          beanplotIndividual(expression_for_pair, neuronal, col2plot, n=n, p=p_val)
        })
      }
      #end_time <- Sys.time()
      #print(end_time - start_time)
      dev.off()
    }
  }

# Removing the testing parts of this function since they are calculated in
# byGRB_list if used so it is added as parameter
beanplotIndividual <- function(expression_for_pair, neuronal, col2plot, n, p) {
  
  par(mar=c(2,1,2,1))
  # start_time <- Sys.time()
  
  
  
  # multiple testing correction, Bonferonni
  p <- n*p

  # ** in p <=0.01, *** if p<=0.001, nothing otherwise
  s <- ifelse(p<=0.001, "***", ifelse(p<=0.01, "**", ifelse(p<=0.05, "*", "")))
  
  expression_for_pair$enh_expression <- log_transform(expression_for_pair$enh_expression, 10)
  expression_for_pair$prom <- log_transform(expression_for_pair$prom, 10)
  
  beanplot(expression_for_pair$prom[expression_for_pair$enh_expression==0], 
           expression_for_pair$prom[expression_for_pair$enh_expression> 0],
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

# Function for calculating empirical p values
empP_diff_in_median_sampling <- function(expression_for_pair, iter=10000) {
  
  # how many tissues have non-zero expression of enhancer
  how_many <- length(which(expression_for_pair$enh_expression>0)) 
  
  # makes iter random samplings of size how_many from enhancer, scores each
  distr <- integer(iter)
  for(i in 1:iter) {
    distr[i] <- median(sample(expression_for_pair$prom, size=how_many, replace=F))
  }
  
  # how many out of iter have lower score than observed -> empirical P
  p_emp <- length(which(distr>=median(expression_for_pair$prom[expression_for_pair$enh_expression>0]))) / iter
  return(p_emp)
}

log_transform <- function(x, base) {
  x <- x+1
  return(log(x, base))
}

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

###
# requires objects mart, byGRBlist, loci.gr and mm10.GRBs.subset
# plots Gviz for each locus and puts it in a .png file
# labels of GRB and locus track can be turned off (others are irresponsice to the flag)



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
