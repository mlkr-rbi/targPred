library(dplyr)
library(readr)
library(stringr)

# take arguments from command line for ppulation
args <- commandArgs(trailingOnly = TRUE)
POP <- args[1]
PLINKDIR <- args[2]
output_file <- args[3]


readCSV <- function(filename, separator) {
    tryCatch({
        df <- read_delim(filename, delim = separator, col_types = cols())
        return(df)
    }, error = function(e) {
        if (grepl("cannot open", e$message)) {
            message(sprintf("ERROR: the file %s does not exist!", filename), appendLF = TRUE)
            quit(status = 1)
        } else if (grepl("no lines available", e$message)) {
            message(sprintf("ERROR: the file %s is empty!", filename), appendLF = TRUE)
            quit(status = 1)
        } else {
            stop(e)
        }
    })
}

frames_list <- list()

print(paste("Processing population:", POP))

ld_files <- list.files(file.path(PLINKDIR, "ld_blocks"), pattern = paste0("^", POP, "_phase3_chr(\\d+|X)\\.ld$"), full.names = TRUE)
print(ld_files)

# 22 plus X chromosome
for (i in c(1:23)) {
    plink_df <- read.table(ld_files[i], header = TRUE)
    plink_grouped_min <- plink_df %>%
        group_by(CHR_A, SNP_A) %>%
        summarise(BP_B_start = min(BP_B), .groups = "drop")
    plink_grouped_max <- plink_df %>%
        group_by(CHR_A, SNP_A) %>%
        summarise(BP_B_end = max(BP_B), .groups = "drop")
    plink_grouped <- left_join(plink_grouped_min, plink_grouped_max, by = c("CHR_A", "SNP_A")) %>%
        mutate(
            CHR = paste0("chr", CHR_A),
            SNPS = as.character(SNP_A)
        ) %>%
        select(CHR, BP_B_start, BP_B_end, SNPS)
    # TODEL
    #write.table(plink_grouped, file.path = , col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    frames_list[[i]] <- plink_grouped
}

plink_allChr <- bind_rows(frames_list)
write.table(plink_allChr, output_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")