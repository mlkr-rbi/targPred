# GWAS table with LD block

import pandas as pd
import numpy as np
#from plydata import define, select
#from plydata.cat_tools import f
#from plydata.pipeline import pipe
import re

# Read in gwas_snps_LDblocks and gwas_associations
gwas_snps_LDblocks = pd.read_csv("out/EUR_phase3_chrALL.bed",sep="\t")
gwas_snps_LDblocks.set_axis(["chromosome","min_loc","max_loc","SNPS"],axis='columns',inplace=True)
gwas_associations = pd.read_csv("/storage/PROJECTS/targPred/dalibor/targPred/plink_jobs/test_plink2/hg37/SNP_2021/gwas_associations.csv")

# Removing "_grch37" from the query snp name
#gwas_snps_LDblocks["SNPS"] = (
#    gwas_snps_LDblocks["query_snp"].str.split(pat = "_", expand = True)[0]
#)

# Creating a column for LDblock that can be easily converted to gr
# Format: chr1:1234-1456
gwas_snps_LDblocks["LDblock"] = (
    gwas_snps_LDblocks["chromosome"] + ":" +
    gwas_snps_LDblocks["min_loc"].astype(str) + "-" +
    gwas_snps_LDblocks["max_loc"].astype(str)
)

# Selecting relevant columns
gwas_associations_LDblocks = (
    gwas_snps_LDblocks[["SNPS", "LDblock"]]
    .merge(gwas_associations, on = "SNPS", how = "left")
)

# Save result to RDS file
gwas_associations_LDblocks.to_csv("gwas_associations_LDblock.csv")
