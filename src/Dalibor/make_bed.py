import os
import sys
import time
import numpy as np
import math
import pandas as pd

def readCSV(filename,separator):
    try:
        np.set_printoptions(threshold=np.inf)

        # Read in DataSet as DataFrame
        df = pd.read_csv(filename,sep=separator, header=0)
    except FileNotFoundError:
        sys.stderr.write(f"ERROR: the file {filename} does not exist!")
        sys.exit()
    except pd.errors.EmptyDataError:
        sys.stderr.write(f"ERROR: the file {filename} is empty!")
        sys.exit()
    return df

for i in range(1,23):
    plink_df = readCSV(f"EUR_phase3_chr{i}.ld","\s+")
    plink_grouped_min = plink_df.groupby(["CHR_A","SNP_A"], as_index=False)["BP_B"].min()
    plink_grouped_max = plink_df.groupby(["CHR_A","SNP_A"], as_index=False)["BP_B"].max()
    plink_grouped = plink_grouped_min.join(plink_grouped_max, lsuffix='_start', rsuffix='_end')
    plink_grouped = plink_grouped.drop(["CHR_A_end","SNP_A_end"],axis=1).astype(str)
    plink_grouped = plink_grouped.rename(columns={"CHR_A_start": "CHR_A", "SNP_A_start": "SNPS"})
    plink_grouped["CHR"] = "chr" + plink_grouped["CHR_A"]
    plink_grouped = plink_grouped[["CHR","BP_B_start","BP_B_end","SNPS"]]
    plink_grouped.to_csv(f"EUR_phase3_chr{i}.bed",index=False,header=False,sep="\t")
