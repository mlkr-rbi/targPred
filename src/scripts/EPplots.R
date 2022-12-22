library("stringr")
library("dplyr")
library("GenomicRanges")
library("RColorBrewer")
library("beanplot")


# Function 1
plotEnhPromBeanPerTrait <- 
  function(byGRBlist,
           i, 
           neuronal=F, 
           sign="permutation",
           loci.gr,
           dir,
           file_suffix="_gene_expression_beanplot_emp.png",
           enhancer_tpm.gr,
           promoter_tpm.df) {
    
    # description: add description and paramater definitions
    
    # Vector of enhancers in a GRB
    enh_vec <- rownames(byGRBlist[[i]][[3]])
    # Generating granges of enhancers in order to find overlaps with trait loci
    if(length(enh_vec) != 0){
      print("No enhancer peaks found.")
      
      enh_vec.gr <- grangesFromDataframe(as.data.frame(enh_vec), colName = enh_vec)
      # Vector of enhancers that overlap trait associated loci
      enh_vec.gr <- subsetByOverlaps(enh_vec.gr, loci.gr)
      
      enh_vec <- mcols(enh_vec.gr)$enh_vec
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

# Function 2
beanplotIndividual <- function(expression_for_pair, neuronal, col2plot, n, p) {
  
  par(mar=c(2,1,2,1))
  # start_time <- Sys.time()
  
  
  
  # multiple testing correction, Bonferonni
  p <- n*p

  s <- ifelse(p<=0.001, "***", ifelse(p<=0.01, "**", ifelse(p<=0.05, "*", "")))
  
  expression_for_pair$enh_expression <- log_transform(expression_for_pair$enh_expression, 10)
  expression_for_pair$prom <- log_transform(expression_for_pair$prom, 10)
  
  beanplot(expression_for_pair$prom[expression_for_pair$enh_expression==0], 
           expression_for_pair$prom[expression_for_pair$enh_expression> 0],
           show.names=F, bw="nrd0", 
           log="", side="both",
           overallline = "median",
           beanlines="median",
           col= list(brewer.pal(8, "Dark2")[8], c(col2plot, "white")),
           what=c(0,1,1,0), border=NA, ylim=c(0,3)) 
  
  legend(x=0.75, y=3.75, legend=s, bty="n", cex=6, text.col="black")
}

# Function 3
log_transform <- function(x, base) {
  x <- x+1
  return(log(x, base))
}

# Function 4
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

# Main 

# Information on loci and traits
loci_in_GRBs.df <- readRDS("data/loci_in_GRBs.df.rds")
loci_in_GRBs.gr <- readRDS("data/loci_in_GRBs.gr.rds")
# byGRB_list
byGRB_list <- readRDS("data/byGRB_list.rds")
# Enhancer and promoter TPMs 
fantom5.enhancers.tpm.gr <- readRDS("data/fantom5.enhancers.tpm.pooled.gr.rds")
promoter_cage_tpm_df <- readRDS("data/promoter_cage_tpm_df_fixed_pooled.rds")

## Creating directory structure
ep_pair_dir <- "data/expression_pairs"

if(!file.exists(ep_pair_dir)){
  dir.create("data/expression_pairs")
}

for (trait in unique(loci_in_GRBs.df$MAPPED_TRAIT)[19:2032]) {
  
  # Generating folder in ep_pair_dir for trait
  trait_dir <- str_replace_all(trait, " |/", "_")
  if(!file.exists(file.path(ep_pair_dir,trait_dir))){
    dir.create(file.path(ep_pair_dir,trait_dir))
  }
  
  # Grbs that are relevant for the trait
  trait_grbs <- 
    filter(loci_in_GRBs.df, MAPPED_TRAIT == trait) %>% 
    pull(grb.coords) %>% 
    unique() %>%
    str_replace_all(pattern = ":|-", replacement = "_") 
  
  for(i_grb in trait_grbs){
    # For continuing grbs
    #for(i_grb in trait_grbs[match(i_grb, trait_grbs):length(trait_grbs)]){
    
    # loci_in_GRB for working trait
    i_loci.gr <- loci_in_GRBs.gr[(mcols(loci_in_GRBs.gr)$MAPPED_TRAIT == trait)]
    plotEnhPromBeanPerTrait(byGRB_list, i=i_grb,
                            # loci associated with trait
                            loci.gr = i_loci.gr,
                            neuronal=F, sign="permutation", 
                            dir=file.path(ep_pair_dir,trait_dir),
                            file_suffix= paste0(trait_dir,"_", i_grb, "_gene_expression_beanplot.png"),
                            enhancer_tpm.gr = fantom5.enhancers.tpm.gr,
                            promoter_tpm.df = promoter_cage_tpm_df)
    
  }
}