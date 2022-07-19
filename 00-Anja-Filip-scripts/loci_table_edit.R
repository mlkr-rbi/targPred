
library(tidyverse)
setwd("/home/filipu/Documents/IRB/GRB/Testing/reactable_testing_new")


grb_dt <- readRDS("02-loci_in_GRB_df-modified_MAPPED_TRAIT.rds")

orig_table <- readRDS("granges_partial_MAPPED_TRAIT_sep.rds")
orig_table_one <- orig_table
orig_table_full <- readRDS("granges_full_MAPPED_TRAIT_sep.rds")


#rename columns 
grb_dt <- rename(grb_dt, Chromosome = chr, SNP = query_snp, Locus_position = "LDblock.0.6", 
                 Study_accession = "STUDY.ACCESSION",  SNP_position = "SNP.position", 
                 Strongest_SNP_risk_allele = "STRONGEST.SNP.RISK.ALLELE", P_value = "P.VALUE", 
                 Mapped_trait = "MAPPED_TRAIT", Reported_genes = "REPORTED.GENE.S.")
orig_table <-  rename(orig_table, Chromosome = chromosome, SNP = query_snp , Locus_position = "LDblock.0.6", 
                      Study_accession = "STUDY.ACCESSION",  SNP_position = "CHR_POS", 
                      Strongest_SNP_risk_allele = "STRONGEST.SNP.RISK.ALLELE", P_value = "P.VALUE", 
                      Mapped_trait = "MAPPED_TRAIT", Reported_genes = "REPORTED.GENE.S.")
grb_dt$Chromosome <- as.character(grb_dt$Chromosome)
grb_dt$SNP_position<- as.character(grb_dt$SNP_position)

#remove NR, NA, intergenic from reported genes and mapped genes, concatenate reported and mapped genes
orig_table <- mutate(orig_table, Reported_genes = gsub("NR", "", Reported_genes))
orig_table <- mutate(orig_table, Reported_genes = gsub("intergenic" | "Intergenic", "", Reported_genes))
orig_table$Reported_genes[is.na(orig_table$Reported_genes)] <- ""
orig_table$MAPPED_GENE[is.na(orig_table$MAPPED_GENE)] <- ""

orig_table$gwas_genes <-  paste(orig_table$MAPPED_GENE, orig_table$Reported_genes, sep= ", ")
#orig_table$gwas_genes <- as.character(orig_table$gwas_genes)
#orig_table$gwas_genes <- data.frame(gwas_genes = vapply(strsplit(orig_table$gwas_genes, ",+"), function(x) paste(unique(x), collapse = ","), character(1)))

#remove columns that are not needed 
grb_dt <- select(grb_dt, -width, -start, -end, -STUDY, -"DISEASE.TRAIT", -REGION, -Reported_genes,
                 -MAPPED_GENE, -SNP_GENE_IDS)
orig_table <- select(orig_table, -STUDY, -"DISEASE.TRAIT", -REGION, -Reported_genes, 
                     -MAPPED_GENE, -SNP_GENE_IDS)



#convert all empty cells in GRBtargets to "unknown", remove "chr" from chromosome number  
grb_dt <- mutate(grb_dt, GRBtargets = gsub('^$', 'unknown', GRBtargets))
grb_dt <- mutate(grb_dt, Chromosome = str_replace(Chromosome, "chr", ""))
#lowercase for all strings in Mapped_trait
grb_dt$Mapped_trait <- tolower(grb_dt$Mapped_trait)
#rename columns
colnames(orig_table)


#remove duplicate rows... note: origin of duplications unknown, one possible non-error duplication is 
#the protein measurement trait where the specific protein association is encoded in a different column
dupl_check_orig <- orig_table 
dupl_check_orig$duplication <- duplicated(dupl_check_orig)
dupl_check_orig_T <- filter(dupl_check_orig, duplication == T)
orig_table <- filter(dupl_check_orig, duplication == F)

dupl_check_dt <- grb_dt 
dupl_check_dt$duplication <- duplicated(dupl_check_dt) 
dupl_check_dt_T <- filter(dupl_check_dt, duplication == T)
grb_dt <- filter(dupl_check_dt, duplication == F)

#from orig_table extract only those rows which are not present in grb_dt
missing_rows_table <- anti_join(orig_table, grb_dt, by = "SNP")

#merge rows extracted in previous line with grb_dt
complete_table <- full_join(grb_dt, missing_rows_table)

#when the GRBtargets cell is NA, put false in inGRB column
complete_table <- mutate(complete_table, inGRB = !is.na(complete_table$GRBtargets))
#remove duplication column
complete_table <- select(complete_table, -duplication)

#change column order for complete table
col_order_comp <- c("Chromosome", "SNP_position", "SNP", "Locus_position","Strongest_SNP_risk_allele", "P_value", "Mapped_trait", "Study_accession",
               "gwas_genes", "GRBtargets", "inGRB")
complete_table <- complete_table[, col_order_comp]

#concatenate Mapped_trait and Locus_position columns to serve as intermediary identifiers 
loci_table <- complete_table
loci_table$merged <- paste(loci_table$Mapped_trait, loci_table$Locus_position, sep= "_")

#change column order for complete table
col_order <- c("merged", "Mapped_trait", "Locus_position", "Chromosome","SNP", "Study_accession", "SNP_position", "Strongest_SNP_risk_allele",
               "P_value", "GRBtargets", "gwas_genes", "inGRB")
loci_table <- loci_table[, col_order]

#add new column "duplication" and mark ALL duplicated instances in "merged" column as TRUE  
loci_table$duplication <- duplicated(loci_table[, 1]) | duplicated(loci_table[, 1], fromLast=TRUE)
#extract only TRUE (duplicated) rows for further inspection 
loci_table_T <- filter(loci_table, duplication == T)

#define function that concatenates, keeps only one instance if duplicated
concat_unique <- function(x){paste(unique(x),  collapse=',')}
#concatenate rows that are "merge" column duplicates, use concat_unique to keep only unique
loci_table_2 <- loci_table %>% 
  group_by(merged) %>% 
  summarise_all(concat_unique)

loci_table_2$gwas_genes <- as.character(loci_table_2$gwas_genes)
#remove duplicates from "gwas_genes" using function below to perform this on strings
loci_table_2$gwas_genes <- data.frame(gwas_genes = vapply(strsplit(loci_table_2$gwas_genes, ",+"), function(x) paste(unique(x), collapse = ", "), character(1)))

#remove unnecessary columns
all_loci <- select(loci_table_2, Chromosome, Locus_position, gwas_genes, Mapped_trait, inGRB, GRBtargets)
all_loci <- as.data.frame(all_loci)

#tidying the cells in gwas_genes
all_loci <- mutate(all_loci, gwas_genes2 = gsub("^,", "", all_loci$gwas_genes2)) 
all_loci$gwas_genes2 <- gsub("^,", "", all_loci$gwas_genes2)

length(unique(complete_table$Mapped_trait))
length(unique(all_loci$Mapped_trait))

saveRDS(complete_table, file="grb_dt_v4.rds")
saveRDS(all_loci, file = "all_loci.rds")
