library(LDlinkR)

wd <- "/home/kprskalo/targPred/LDproxy"


dir <- file.path(wd,"LDproxy")
dir.create(dir)
setwd(dir)

failed_snps <- readRDS(file.path("LDproxy_failed_snps.rds"))

ldlink_token <- "be94d998f060"

for(rs in failed_snps){
  
  # Running LDproxy_batch 
  # tryCatch skips errors that occur because of snps missing from LDproxy databases
  tryCatch(LDproxy_batch(rs, pop = "CEU", r2d = "r2", token = ldlink_token, append = F),
           error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}