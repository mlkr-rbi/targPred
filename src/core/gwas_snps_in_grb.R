
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)
library(httr)
library(jsonlite)

core_dir <- "src/core" 

source(file.path(core_dir, "utils.R"))

path_file <- "src/core/paths.txt"
paths <- set_paths(path_file)
core_dir <- paths$core_dir
work_dir <- paths$work_dir
data_dir <- paths$data_dir
empP_dir <- paths$empP_dir
gwas_year <- paths$gwas_year
LDblock_path <- paths$LDblock_path
byGRBlist_name <- paths$byGRBlist_name




setwd(work_dir)

### Script inputs ###

## GRBs with target genes
grbs <- readRDS(file.path(data_dir, "processed", "GRBs_with_targets.gr.rds"))

## Load EFO mappings from GWAS - this is downloaded from the FTP site as same as associations file
efo_mappings <- read.csv(file.path(data_dir, "inputs", paste0("GWAS_",gwas_year), "gwas-efo-trait-mappings.tsv"), sep = "\t", header = TRUE)

## GWAS association file
gwas_path <- file.path("data/processed", "gwas_associations_2026.csv")


### Script outputs ###

## Log file
log_fn <- file.path(data_dir, "processed/logs", "gwas_snps_in_grb.txt")

## GWAS associations with LD blocks and GRB targets
gwas_fn <- file.path(data_dir,"outputs", "gwas_associations.csv")

## Bed file with loci in GRBs for visualizations in genome browser
loci_in_GRBs_bed_fn <- file.path(data_dir, "outputs", "loci_in_GRBs.bed")

## Dataframe with loci in GRBs for further analyses
loci_in_GRBs_df_fn <- file.path(data_dir, "processed","loci_in_GRBs.df.csv")

## Granges object with SNPs in GRBs for further analyses
 snps_in_GRBs_gr_fn <- file.path(data_dir, "processed","snps_in_GRBs.gr.rds")

## Mart result for SNPs in GWAS
mart_result_fn <- file.path(data_dir,"inputs", "GWAS_SNPs_mart_result.csv")

# Web table with all the information for the web app
web_table_fn <- file.path(data_dir,"outputs", "web_table.csv")

# SQL database with all the information for the web app
sql_db_fn <- file.path(data_dir,"outputs", paste0("targpred", gwas_year, ".db"))


gwas_associations <-
  read.csv(file = gwas_path, header = TRUE, stringsAsFactors = FALSE, row.names = NULL) %>%
  unique()

# Cat the dimensions of gwas_associations to log
cat("Dimensions of gwas_associations after EUR filtering:", dim(gwas_associations), "\n", file = log_fn, append = TRUE)

gwas_snps_LDblocks <- read.csv(file.path(LDblock_path),sep="\t",header=FALSE)
colnames(gwas_snps_LDblocks) <- c("chromosome","min_loc","max_loc","SNPS")
gwas_snps_LDblocks <-
  gwas_snps_LDblocks %>%
  # Creating a column for LDblock that can be easily converted to gr 
  # Convert chr23 to chrX
    mutate(chromosome = case_when(chromosome == "chr23" ~ "chrX",
                                    chromosome == "chr24" ~ "chrY",
                                    TRUE ~ as.character(chromosome))) %>%
  # Format: chr1:1234-1456
  mutate(LDblock = str_c(chromosome,":",min_loc,"-",max_loc)) %>%
  dplyr::select(SNPS, LDblock)

# Save the number of LD blocks to log
cat("Number of LD blocks from GWAS EUR:", nrow(gwas_snps_LDblocks), "\n", file = file.path(data_dir, "processed/logs", "gwas_snps_in_grb.txt"), append = TRUE)



# Join LD blocks with GWAS associations
gwas <- left_join(gwas_snps_LDblocks, gwas_associations, by="SNPS") %>%
    distinct(.keep_all = TRUE)

# Cat the dimensions of gwas after joining with LD blocks to log
cat("Dimensions of gwas after joining with LD blocks:", dim(gwas), "\n", file = log_fn, append = TRUE)


 
### Modifying genes from gwas data #####
# Gwas catalog stores the putative target gene of an SNP in two columns: (from https://www.ebi.ac.uk/gwas/docs/fileheaders ) REPORTED GENE(S): Gene(s) reported by author MAPPED GENE(S): Gene(s) mapped to the strongest SNP. If the SNP is located within a gene, that gene is listed. If the SNP is located within multiple genes, these genes are listed separated by commas. If the SNP is intergenic, the upstream and downstream genes are listed, separated by a hyphen. - Combine the two columns in a single column GWAS_genes

non_gene_names <- c("",NA,"intergenic","Intergenic", "NR")

gwas <- gwas %>%
    # Creating a temporary column that contains all genes comma sep
    mutate(temp = case_when(is.element(REPORTED.GENE.S., non_gene_names) ~ MAPPED_GENE,
                                  is.element(MAPPED_GENE, non_gene_names) ~ REPORTED.GENE.S.,
                                  TRUE ~ str_c(REPORTED.GENE.S., MAPPED_GENE, sep = ", "))) %>%
    # Removing repeated genes (eg. SKIV2L, SKIV2L, SKIV2L, SKIV2L, SKIV2L, SKIV2L)
    mutate(temp1 = map_chr(str_split(temp, ", |;| - |-"), ~str_c(unique(.x),collapse = ", "))) %>%
    mutate(GWAS_GENES = map_chr(str_split(temp1, ", "), ~str_c(.x[!.x == "NA"],collapse = ", "))) %>%
    dplyr::select(-c(temp, temp1, REPORTED.GENE.S., MAPPED_GENE)) 


### Separating Traits ####
# Each investigated SNP has an associated mapped trait, or multiple. Since our goal is to explore the relationship between a SNP and a single trait row separation based on trait is required. Trait names are stored in MAPPED_TRAIT and the URI that contain unique disease/trait ontology id in MAPPED_TRAIT_URI. Traits are comma separated but some trait names contain commas. To avoid error prone case by case exceptions, trait separation will be based on the URI and trait names will be cross mapped from the ontology.
efo_mappings_sub <- efo_mappings %>%
    select(EFO.term, EFO.URI) %>%
    unique()

# Separate GWAS rows by trait using EFO mappings
gwas <- gwas %>%
    tidyr::separate_rows(MAPPED_TRAIT_URI, sep=",") %>%
    mutate(MAPPED_TRAIT_URI = str_trim(MAPPED_TRAIT_URI)) %>%
    left_join(efo_mappings_sub, by = c("MAPPED_TRAIT_URI" = "EFO.URI")) %>%
    select(-MAPPED_TRAIT) %>%
    select(everything(), MAPPED_TRAIT = EFO.term)
     # Remove rows with NA in MAPPED_TRAIT_URI or MAPPED_TRAIT


# Cat the dimensions of gwas after joining with LD blocks and separating traits
cat("Dimensions of gwas after separating traits:", dim(gwas), "\n", file = file.path(data_dir, "processed/logs", "gwas_snps_in_grb.txt"), append = TRUE)


# Drop some columns that are not needed for the web table and further analyses
gwas <- gwas %>%
    dplyr::select(-c(STUDY, DISEASE.TRAIT, SNP_GENE_IDS, STRONGEST.SNP.RISK.ALLELE, MAPPED_TRAIT_URI, OR.or.BETA, X95..CI..TEXT.)) %>%
    unique() 


cat("Dimensions of gwas after dropping unnecessary columns and applying unique:", dim(gwas), "\n", file = log_fn, append = TRUE)
  
write.csv(gwas, file = gwas_fn, row.names = FALSE)


#### Loci in GRBs ####

## gwas_gr stores SNP LD blocks in GenomicRanges object
gwas.gr <- grangesFromDataframe(gwas, LDblock)

# Overlaping gwas loci with GRBs
loci_in_GRBs.gr <- subsetByOverlaps(gwas.gr, grbs)
ovlp <- findOverlaps(loci_in_GRBs.gr, grbs)

# Number of loci in GRBs
cat("Number of loci in GRBs:", length(unique(queryHits(ovlp))), "\n", file = log_fn, append = TRUE)

# Number of unique GRBs that overlap with loci out of total GRBs
cat("Number of unique GRBs that overlap with loci:", length(unique(subjectHits(ovlp))), "out of", length(grbs), "\n", file = log_fn, append = TRUE)

# Adding GRB coordinates to each locus as character
loci_in_GRBs.gr$grb.coords[queryHits(ovlp)] <-
    as.character(grbs[subjectHits(ovlp)])

# Adding GRB target genes 
loci_in_GRBs.gr$GRB_GENES[queryHits(ovlp)] <- 
    grbs[subjectHits(ovlp)]$geneName

# Saving loci in GRB as bed file for visualizations in Genome browser
export.bed(loci_in_GRBs.gr, con = loci_in_GRBs_bed_fn)

loci_in_GRBs.df <- mcols(loci_in_GRBs.gr) %>% as.data.frame()


write.csv(loci_in_GRBs.df, file = loci_in_GRBs_df_fn, row.names = FALSE )

#LOG: number of unique SNPs in GWAS
cat("Number of unique SNPs in GWAS:", length(unique(gwas$SNPS)), "\n", file = log_fn, append = TRUE)
#LOG: number of unique SNPs in GRBs
cat("Number of unique SNPs in GRBs:", length(unique(loci_in_GRBs.df$SNPS)), "\n", file = log_fn, append = TRUE)

# SNPs locations 
if (file.exists(mart_result_fn)) {
    mart_result <- read.csv(file = mart_result_fn)
} else {
    getBM_retry <- function(...) {
  for (i in 1:5) {
    tryCatch({
      return(getBM(...))
    }, error = function(e) {
      message(sprintf("Attempt %d failed: %s", i, e$message))
      Sys.sleep(10)
    })
  }
  stop("All attempts to query biomaRt failed.")
  }
    ensembl <- useMart(biomart = "ENSEMBL_MART_SNP", host = "https://grch37.ensembl.org", dataset="hsapiens_snp")
    # Split rs_numbers into batches of 5000
    batch_size <- 5000
    rs_numbers <- unique(gwas$SNPS)
    mart_result <- data.frame()
    for (i in seq(1, length(rs_numbers), by = batch_size)) {
    message(sprintf("Processing batch %d to %d of %d", i, min(i+batch_size-1, length(rs_numbers)), length(rs_numbers)))
    batch <- rs_numbers[i:min(i+batch_size-1, length(rs_numbers))]
    res <- getBM_retry(
        attributes = c('refsnp_id', 'chr_name', 'chrom_start'),
        filters = 'snp_filter',
        values = batch,
        mart = ensembl
    )
    mart_result <- rbind(mart_result, res)
    # Exclude non-cannonical chromosomes (eg. chrUn, chrM, chrX_random, etc.)
    mart_result <- mart_result %>%
    filter(!str_detect(chr_name, pattern = "^H"))

    # Save the mart_result to a csv file for future use
    write.csv(mart_result, file = mart_result_fn, row.names = FALSE)
    }
    
    # Check which SNPs in gwas$SNPS are not present in mart_result - those are likely archived SNPs that are not in dbSNP anymore, or SNPs that do not have a position in the database
    missing_snps <- setdiff(gwas$SNPS, mart_result$refsnp_id)

    # Use httr and jsonlite to query dbSNP REST API for the missing SNPs to get their positions
    if (length(missing_snps) > 0) {
        message(sprintf("Querying dbSNP REST API for %d missing SNPs", length(missing_snps)))
        missing_snps_df <- data.frame(refsnp_id = missing_snps, cannonical_id = NA, chr_name = NA, chrom_start = NA)
        for (i in 1:nrow(missing_snps_df)) {
            # print progress every 10 SNPs
            if (i %% 10 == 0) {
                message(sprintf("Processing missing SNP %d of %d", i, nrow(missing_snps_df)))
            }
            snp_id <- missing_snps_df$refsnp_id[i]
            url <- paste0("https://rest.ensembl.org/variant_recoder/human/", snp_id, "?content-type=application/json")
            response <- httr::GET(url)
            data <- fromJSON(rawToChar(response$content))
            # Get the cannonical id
            missing_snps_df$cannonical_id[i] <- unlist(data)["TT.id"]
        }

        # save missing_snps_df to a csv file for future use
        write.csv(missing_snps_df, file = file.path(data_dir,"processed", "missing_snps_df.csv"), row.names = FALSE)
        # For the cannonical ids, query biomaRt again to get their positions
        res <- getBM(
            attributes = c('refsnp_id', 'chr_name', 'chrom_start'),
            filters = 'snp_filter',
            values = missing_snps_df$cannonical_id,
            mart = ensembl
        )
        # Combine the results with the original mart_result
        mart_result <- rbind(mart_result, res)
        write.csv(mart_result, file = mart_result_fn, row.names = FALSE)

    }
}

#LOG: how many of the SNPs in gwas$SNPS and loci_in_GRBs.df$SNPS are not present in mart_result
cat("Number of SNPs in gwas$SNPS that are not present in mart_result:", length(setdiff(gwas$SNPS, mart_result$refsnp_id)), "\n", file = log_fn, append = TRUE)
cat("Number of SNPs in loci_in_GRBs.df$SNPS that are not present in mart_result:", length(setdiff(loci_in_GRBs.df$SNPS, mart_result$refsnp_id)), "\n", file = log_fn, append = TRUE)


# Extracting SNPs that are in GRBs and creating a GRanges object for them
snps_in_GRBs.gr <- filter(mart_result, refsnp_id %in%      loci_in_GRBs.df$SNPS) %>%
# Rename columns
    select(SNPS = refsnp_id, CHR = chr_name, POS = chrom_start) %>%
    # Creating a GRanges object
    makeGRangesFromDataFrame(seqnames.field = "CHR", start.field = "POS", end.field = "POS", keep.extra.columns = TRUE)

# Concetinating  to the final object
## Changing type of seqnames from 1,2,3,... to chr1, chr2, chr3,...
 seqlevelsStyle(snps_in_GRBs.gr) <- "UCSC"
# 
# # ?Maybe combine with loci_in_GRBs.df for extra metadata
saveRDS(snps_in_GRBs.gr, file = snps_in_GRBs_gr_fn)


### Prepare web_table ###
web_table <- gwas %>%
    left_join(select(loci_in_GRBs.df, SNPS, grb.coords, GRB_GENES), by = c("SNPS")) %>%
    distinct(.keep_all = FALSE) 

web_table <- web_table %>%
    # Rename columns for web table
    dplyr::select(
        Chromosome = CHR_ID,
        SNP = SNPS,
        Locus_position = LDblock,
        P_value = P.VALUE,
        Mapped_trait = MAPPED_TRAIT,
        Study_accession = STUDY.ACCESSION,
        GRB_position = grb.coords,
        GRB_targets = GRB_GENES,
        GWAS_genes = GWAS_GENES
    ) %>%
    # Add a column indicating whether the SNP is in a GRB or not
    mutate(inGRB = !is.na(GRB_position)) %>%
    # If the SNP is in a GRB but has no known targets, set GR_targets to "unknown"
    mutate(GRB_targets = if_else((GRB_targets == "NA" & inGRB == TRUE), "unknown", GRB_targets))

# Add SNP positions to the web table by joining with mart_result
web_table <- web_table %>%
    left_join( select(mart_result, refsnp_id, chrom_start), by = c("SNP" = "refsnp_id")) 

web_table <- web_table %>%
# Rename chrom_start to SNP_position and reorder columns using select
    select("Chromosome","SNP","Locus_position", SNP_position = "chrom_start","P_value","Mapped_trait","Study_accession","GRB_position","GRB_targets","GWAS_genes","inGRB")


write.csv(web_table, file = web_table_fn, row.names = FALSE)

conn <- dbConnect(RSQLite::SQLite(), sql_db_fn)
dbWriteTable(conn, "grb_table", web_table, overwrite = TRUE)
dbWriteTable(conn, "loci_table", loci_in_GRBs.df, overwrite = TRUE)
dbDisconnect(conn)

#LOG number of rows in web_table
cat("Number of rows in web_table:", nrow(web_table), "\n", file = log_fn, append = TRUE)
#LOG number of unique SNPs in web_table
cat("Number of unique SNPs in web_table:", length(unique(web_table$SNP)), "\n", file = log_fn, append = TRUE)
#LOG number of unique traits in web_table
cat("Number of unique traits in web_table:", length(unique(web_table$Mapped_trait)), "\n", file = log_fn, append = TRUE)
#LOG number of unique loci in web_table
cat("Number of unique loci in web_table:", length(unique(web_table$Locus_position)), "\n", file = log_fn, append = TRUE)


