library(LDlinkR)

wd <- getwd()

gwas_snps_ls <- readRDS(file.path(wd,"gwas_snps_ls.rds"))

ldlink_token <- "be94d998f060"


for(chr in names(gwas_snps_ls)[]){
  

  # Creating a dir for each chromosome
  dir <- file.path(wd,chr)
  dir.create(dir)
  setwd(dir)
  
  # Nested for loop for skipping errors
  for(rs in gwas_snps_ls[[chr]]$SNPS){
    
  # Running LDproxy_batch 
  tryCatch(LDproxy_batch(rs, pop = "CEU", r2d = "r2", token = ldlink_token, append = F),
           error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
}