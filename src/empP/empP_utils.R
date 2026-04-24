library(GenomicRanges)
library(dplyr)
library(stringr)
library(stats)


calcWilcox <- function(row, y) {
  tryCatch(
    {
      # Apply your wilcox.test or other function here
      result <- wilcox.test(row, y)$statistic
    },
    error = function(e) {
      # Return a default value if an error occurs
      result <- NA
    }
  )
  return(result)
}

calcTTest <- function(row, y) {
  tryCatch(
    {
      # Apply your wilcox.test or other function here
      result <- t.test(row, y)$statistic
    },
    error = function(e) {
      # Return a default value if an error occurs
      result <- NA
    }
  )
  return(result)
}

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
  
  # message(head(data_r))
  
  # Granges object is needed to detect SNP loci that overlap 2 GRBs
  data_r <-
    makeGRangesFromDataFrame(data_r,
                             seqnames.field ="CHR",
                             start.field = "start",
                             end.field = "end",
                             keep.extra.columns = keep_metadata)
  
  
  return(data_r)
}

# OG - Anjas function
empP_diff_in_median_sampling <- function(expression_for_pair, iter=10000) {
  
  # how many tissues have non-zero expression of enhancer
  how_many <- length(which(expression_for_pair$enh_expression>0)) 
  
  # makes iter random samplings of size how_many from enhancer, scores each
  distr <- integer(iter)
  for(i in 1:iter) {
    distr[i] <- median(sample(expression_for_pair$prom, size=how_many, replace=F))
  }
  
  median_ep <- median(expression_for_pair$prom[expression_for_pair$enh_expression>0])
  # how many out of iter have lower score than observed -> empirical P
  p_emp <- length(which(distr>=median_ep)) / iter


  return(c(how_many, p_emp, median_ep, median(distr)))
}


# Preparing empirical p value
calcEmpP <- function(expression_for_pair, iter=10000) {
  #' Calculate Empirical p-Values for Enhancer Analysis
  #'
  #' This function calculates empirical p-values for enhancer analysis based on random sampling and statistical testing.
  #'
  #' @param expression_for_pair A data frame containing enhancer and promoter expression data.
  #' @param iter Number of random samplings to perform.
  #'
  #' @return A numeric vector containing the number of enhancers with non-zero expression, one-tailed and two-tailed empirical p-values for median, t-test, and Wilcoxon test.
  #'
  #' @details The function calculates empirical p-values by performing random sampling of tissues with active enhancer expression. It then compares these samples to the overall distribution and calculates p-values.
  #'
  #' @examples
  #' # Example usage:
  #' result <- calcEmpP(expression_for_pair, iter = 10000)
  #' 
  
  # how many tissues have non-zero expression of enhancer
  e_active <- length(which(expression_for_pair$enh_expression>0))

  prom_E0 <- expression_for_pair$prom[expression_for_pair$enh_expression==0]
  prom_E1 <- expression_for_pair$prom[expression_for_pair$enh_expression>0]
  
  # makes iter random samplings of size how_many from enhancer, scores each
  distr_median <- integer(iter)
  distr_ttest <- integer(iter)
  # distr_wtest <- integer(iter)

  # for(i in 1:iter) {
  #   # sampling tissues - number is equal to number of tissues with active enhancer

  #   distr <- sample(expression_for_pair$prom, size=e_active, replace=F)
  #   # Calculating median and t-test statistics (z score) for the sample
  #   distr_median[i] <- median(distr)
  #   distr_ttest[i] <- tryCatch({t.test(distr, prom_E0)$statistic}, error = function(e) {
  #     # Return a default value if an error occurs
  #     NA})
  #   distr_wtest[i] <- tryCatch({wilcox.test(distr, prom_E0)$statistic}, error = function(e) {
  #     # Return a default value if an error occurs
  #     NA})
  # }

  # Sample expression_for_pair$prom e_active times without replacement and replicate it iter times
  distr_ls <- replicate(iter, sample(expression_for_pair$prom, size = e_active, replace = FALSE), simplify = TRUE)
  distr_mat <- as.data.frame(t(distr_ls))

  # Calculating median
  distr_median <- apply(distr_mat, 1, median)

  # Calculating ttest
  # distr_ttest <- apply(distr_mat, 1, calcTTest, y = prom_E0)
  
  
  # Calculating median and t-test statistics (z score) for active enhancers
  median_baseline <- median(prom_E1)
  # ttest_baseline <- tryCatch({t.test(prom_E1, prom_E0)$statistic}, error = function(e) {
  #     # Return a default value if an error occurs
  #     NA})
  # wtest_baseline <- tryCatch({wilcox.test(prom_E1, prom_E0)$statistic}, error = function(e) {
  #     # Return a default
  #     NA})
  
  # Calculating p values for onetailed and two tailed
  # Numerator is the number of observations in our tested distribution is higher than sampled
  # Denominator is the number of iterations or replicate samples
  r_median <- length(which(distr_median>=median_baseline))
  p_emp_median_oneTailed <- (r_median+1) / (iter+1)
  r_median <- ifelse(r_median < iter/2, r_median, iter-r_median) 
  p_emp_median_twoTailed <- ((r_median+1) / (iter+1))*2



  # For t_test
  # r_ttest <- length(which(distr_ttest>=ttest_baseline))
  # p_emp_ttest_oneTailed <- (r_ttest+1) / (iter+1)
  # r_ttest <- ifelse(r_ttest < iter/2, r_ttest, iter-r_ttest)
  # p_emp_ttest_twoTailed <- ((r_ttest+1) / (iter+1))*2
  
  # For w_test
  # r_wtest <- length(which(distr_wtest>=wtest_baseline))
  # p_emp_wtest_oneTailed <- (r_wtest+1) / (iter+1)
  # r_wtest <- ifelse(r_wtest < iter/2, r_wtest, iter-r_wtest)
  # p_emp_wtest_twoTailed <- ((r_wtest+1) / (iter+1))*2

  # p_emp_wtest_oneTailed <- NA
  # p_emp_wtest_twoTailed <- NA

  # out_vec <- c(e_active, p_emp_median_oneTailed, p_emp_median_twoTailed, p_emp_ttest_oneTailed, p_emp_ttest_twoTailed, p_emp_wtest_oneTailed, p_emp_wtest_twoTailed)
  # names(out_vec) <- c("Number_of_enhancers", "empP_median_oneTailed", "empP_median_twoTailed", "empP_ttest_oneTailed","empP_ttest_TwoTailed", "empP_wtest_oneTailed", "empP_wtest_TwoTailed")

  # Only median p values are calculated for now, t-test and w-test are commented out
  out_vec <- c(e_active, p_emp_median_oneTailed, p_emp_median_twoTailed)
  names(out_vec) <- c("Number_of_enhancers", "empP_median_oneTailed", "empP_median_twoTailed")

  return(out_vec)
}


# Alternative calcEmpP using median difference permutation test
calcEmpP_alt <- function(expression_for_pair, iter=10000) {
  #' Calculate Empirical p-Values Using Median Difference Permutation Test
  #'
  #' This function calculates empirical p-values for enhancer analysis based on permutation testing of median differences.
  #' It compares the observed difference in medians between tissues with active enhancers vs inactive enhancers
  #' to a null distribution generated by permutation.
  #'
  #' @param expression_for_pair A data frame containing enhancer and promoter expression data with columns: 
  #'   enh_expression (enhancer expression) and prom (promoter expression).
  #' @param iter Number of permutations to perform (default = 10000).
  #'
  #' @return A numeric vector containing the number of tissues with active enhancers, 
  #'   one-tailed and two-tailed empirical p-values based on median difference.
  #'
  #' @details This function performs a permutation test on the difference in medians between
  #'   two groups (active vs inactive enhancer tissues). The observed test statistic is compared
  #'   to the permutation distribution to calculate p-values.
  #'
  #' @examples
  #' # Example usage:
  #' result <- calcEmpP_alt(expression_for_pair, iter = 10000)
  
  # Extract tissues with active and inactive enhancer expression
  prom_E0 <- expression_for_pair$prom[expression_for_pair$enh_expression == 0]
  prom_E1 <- expression_for_pair$prom[expression_for_pair$enh_expression > 0]
  
  e_active <- length(prom_E1)
  
  # Observed test statistic: difference in medians
  obs_stat <- median(prom_E1) - median(prom_E0)
  
  # Combine both groups for permutation
  combined <- c(prom_E1, prom_E0)
  n1 <- length(prom_E1)
  n  <- length(combined)
  
  # Perform permutation test
  perm_stats <- numeric(iter)
  for (b in seq_len(iter)) {
    idx     <- sample.int(n, size = n1, replace = FALSE)
    x_perm  <- combined[idx]
    y_perm  <- combined[-idx]
    perm_stats[b] <- median(x_perm) - median(y_perm)
  }
  
  # Calculate one-tailed p-value (greater than or equal to observed)
  r_one <- sum(perm_stats >= obs_stat)
  p_emp_oneTailed <- (r_one + 1) / (iter + 1)
  
  # Calculate two-tailed p-value (absolute value)
  r_two <- sum(abs(perm_stats) >= abs(obs_stat))
  p_emp_twoTailed <- (r_two + 1) / (iter + 1)
  
  # Return in same format as calcEmpP
  out_vec <- c(e_active, p_emp_oneTailed, p_emp_twoTailed)
  names(out_vec) <- c("Number_of_enhancers", "empP_median_oneTailed", "empP_median_twoTailed")
  
  return(out_vec)
}


EffectSizeOfEnhancerActivityCohens <- function(expression_for_pair){
  #' Calculate Cohen's d Effect Size for Enhancer Activity
  #'
  #' This function calculates the Cohen's d effect size for enhancer activity by comparing two groups based on enhancer expression.
  #'
  #' @param expression_for_pair A data frame containing enhancer and promoter expression data.
  #'
  #' @return The Cohen's d effect size value, which indicates the standardized mean difference between tissues with and without active enhancer expression.
  #'
  #' @details Cohen's d is a measure of the standardized effect size, indicating how far apart the means of two groups are in terms of standard deviations. A larger Cohen's d suggests a larger effect size.
  #'
  #' @examples
  #' # Example usage:
  #' effect_size <- EffectSizeOfEnhancerActivityCohens(expression_for_pair)
  group1 <- expression_for_pair$prom[expression_for_pair$enh_expression>0]
  group2 <- expression_for_pair$prom[expression_for_pair$enh_expression==0]

  # Calculate the mean difference
  mean_diff <- mean(group1) - mean(group2)

  # Calculate the pooled standard deviation
  n1 <- length(group1)
  n2 <- length(group2)
  var1 <- var(group1)
  var2 <- var(group2)
  pooled_sd <- sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

  # Calculate Cohen's d
  cohen_d <- mean_diff / pooled_sd
  return(cohen_d)  # prints the Cohen's d value

}



EnhExpression <- function(fantom5.enh.gr, enh.gr){
  return(apply(as.data.frame(mcols(subsetByOverlaps(fantom5.enh.gr, enh.gr))), 2, sum))
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


CreateEnhPromCombinationsPerChromosome <- function(fantom5.enh.gr, TSS.gr, promoter.df, chr, max_distance = 1000000, metrics_columns, create_new = TRUE, filename = ""){
  #' Create Enhancer-Promoter Combinations Per Chromosome
#'
#' This function generates enhancer-promoter combinations for a specific chromosome based on given genomic data. It calculates the distances between enhancers and promoters, filtering by a maximum distance threshold.
#'
#' @param fantom5.enh.gr A GRanges object representing enhancer locations.
#' @param TSS.gr A GRanges object representing transcription start site (TSS) locations.
#' @param promoter.df A data frame containing promoter information.
#' @param chr The chromosome identifier for which to create combinations.
#' @param max_distance The maximum distance threshold for enhancer-promoter combinations.
#' @param metrics_columns A character vector specifying the names of columns to include in the resulting data.
#' @param create_new Logical, if TRUE (default), create new combinations; if FALSE, read combinations from an existing file.
#' @param filename The filename to read existing combinations from if create_new is set to FALSE.
#'
#' @return A data frame containing enhancer-promoter combinations for the specified chromosome. The columns are specified by the metrics_columns parameter.
#'
#' @details This function creates enhancer-promoter combinations for a specific chromosome by merging enhancer and TSS data. It filters the combinations based on a maximum distance threshold (max_distance). If create_new is set to TRUE, it generates new combinations; otherwise, it reads existing combinations from the provided filename.
#'
#' @examples
#' # Example usage:
#' combinations <- CreateEnhPromCombinationsPerChromosome(fantom5.enh.gr, TSS.gr, promoter.df, chr = "chr1", max_distance = 1000000, metrics_columns = c("enh_pos", "gene", "transcript_position", "EP_dist"))

  if(!create_new){
    empP_tb_res <- read.table(file= file.path(filename), header = TRUE)
    head(empP_tb_res)
    empP_tb_res <- empP_tb_res %>%
      dplyr::select(-chr)

    return(empP_tb_res)
  }

  # Subseting by chromosome
  enh_chr <- fantom5.enh.gr[seqnames(fantom5.enh.gr) == chr]
  transc_chr <- TSS.gr[seqnames(TSS.gr) == chr]

  # Preping some additional info
  transc_chr$gene <- paste0(transc_chr$ensembl_gene_id,":",transc_chr$ensembl_transcript_id)
  transc_chr$chr <- seqnames(transc_chr)
  transc_chr$transcript_position <- start(transc_chr)

  # remove transcir

  # Removing transcripts for which we do not have cage info
  transc_chr <- transc_chr[is.element(transc_chr$gene, rownames(promoter.df))]

  # Converting to dfs for merge
  enh_chr.df <- as.data.frame(granges(enh_chr))

  empP_tb_res <- merge(enh_chr.df, mcols(transc_chr))

  # Filtering by max dist that will be tested
  empP_tb_res <- as.data.frame(empP_tb_res) %>%
    # end - start represnts middle of enhancers
    dplyr::mutate(EP_dist = round(abs((end + start)/2 - transcript_position), digits = 0)) %>%
    dplyr::filter( EP_dist <= max_dist) %>%
    # dplyr:: mutate(enh_pos = paste0(seqnames, ":", start,"-",end,":",strand)) %>%
    dplyr:: mutate(enh_pos = paste0(seqnames, ":", start,"-",end)) %>%
    dplyr::select(enh_pos, gene, transcript_position, EP_dist)

  return(empP_tb_res)

}