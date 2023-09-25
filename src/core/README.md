
# targPred workflow

This workflow includes calculating LD blocks from GWAS data, assigning GWAS and GRB target genes to SNPs and creating plots 


## 1) Preparing 1000Genomse data for calculating LD blocks

Follow workflow on https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html

Download hg37 data from https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3

Be careful to select the correct ensemble (hg19/GRCh27)
Create folder workdir/input and go thither

## 2)	Preprocessing 1000genomes data
```
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen
```

## 3) Select the subpopulation
This selects the correct population (EUR) and generates subpopulations
It also transforms the output to Plink1 format.
To do so run this job
```
qsub plink_job.sge
```
## 4) Calculate the LD blocks
This calculates linkage disequillibrium blocks. It also transforms the output to GRanges.
```
qsub plink1_job.sge
```
# 5) Create bed files from plink output
```
python make_bed.py
```

6) Concatenate files
Concatenate all .bed files using `cat` command to out/EUR_phase3_chrALL.bed

7) Add GWAS information
Run this script to generate gwas_associations_LDblock.csv
```
gwas_associations.py
```

## 8) Overlapping GWAS variants with GRBs
This runs `gwas_snps_in_grb.R` or `gwas_snps_in_grb_no_kurtosis.R`

You need to select kurtosis or no_kurtosis by changing the script name in .sge script.
```
qsub make_gwas_snp_in_grb.sge
```
## 9) Create directories for plots
```
mkdir plots
mkdir plots/expression_pairs
mkdir plots/gviz_plots
```

## 10) Make Gviz plots
These plots represent genomic regions of interest
```
qsub make_gviz.sge
```

## 11) Make beanplots
Beanplots show transcirpt expression depending on corresponding nhancer activity and empirical p value association. 
```
qsub make plots.sge
```
