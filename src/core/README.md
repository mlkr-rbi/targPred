
# targPred workflow

This workflow includes calculating LD blocks from GWAS data, assigning GWAS and GRB target genes to SNPs and creating plots.

The scripts with .sge extensions are to be used with the Sun Grid Engine in a job scheduling HPC cluster setting.

## 1) Create your working directory

```
mkdir my_work_dir
cd my_work_dir
```

## 2) Copy the SGE scripts and paths.txt from the source code directory into your working directory

```
cp /path/to/your/targPred/source/*.sge . 
cp /path/to/your/targPred/source/paths.txt . 
```

After that you need to edit the paths.txt file with the required paths and filnames for targPred to function properly. The file contains keywords and values. The lines containing keywords start with the $ sign, while the following line must contain the appropriate value for that keyword. Any line beginning with a # sign will be interpreted as a comment and ignored.

Meaning of the keywords:

- $core: path to your targPred code directory where the scripts are located
- $data: path to the directory containing data, organized as follows
   1. inputs: directory containing EFO\_terms.rds, GWAS catalog, hg19.canFam3.98.50.processed.threshold.grbs.bed  and hg19.canFam3.tsv 
   2. processed: directory containing the directory with the empP files, genes.hg19.ensembl.rds, TSSs.perTranscript.gr.rds, promoter\_cage\_tpm\_df\_fixed\_pooled.rds and fantom5.enhancers.tpm.pooled.gr.rds. Some intermediary files are also stored here.
   3. outputs: the resulting files, such as byGRB\_list.rds are stored here.
- $work: path to your working directory, where you have copied the SGE scripts and paths.txt.
- $empP\_dir: directory name with your empP files
- $empP: file prefix of empP files
- $year: year of the GWAS catalog edition
- $LDblock\_path: full path leading to the bed file with LD blocks, resulting form the PLINK analysis - the name of the file is EUR_LD_blocks_allChr.bed
- $byGRB\_name: filename to save beGRBi\_list 

The SGE scripts can be adapted for usage on other job scheduling systems such as Slurm or PBS. 


## Downloading GWAS data and calculating LD blocks using 1000Genomes data



To process the GWAS catalog, use the `process_gwas_catalog.R` script. This script downloads the latest GWAS catalog (e.g., the `/2026/02/17/` release), parses the relevant SNP and trait information, and formats the data for downstream analysis. 


The script `run_PLINK.sge` follows the workflow, as outlined in `https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html`.l It downloads the hg19 data from `https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3`, selects the desired population, transforms the data into PLINK1 format and calculates the linkage disequilibrium blocks. 

The user should add the absolute or relative paths to directories containing PLINK1, PLINK2 and Python executables, as well as the path to the text file containing the list of SNPs for the LD block calculation. The user can change the PLINK settings for LD block calculations directly in the script.


The script can afterwards be submitted to the queueing system through the following command:

```
squeue run_PLINK.sh
```

The bed file with the linkage disequilibrium blocks can be found in the newly created directory named ld\_blocks. The script produces a file with LD blocks for each chromosome separately, as well as their concatenated version named "${POP}.bed", where ${POP} denotes the desired population. Bed file with all combined chromosomes is stored in data/processed directory: EUR_LD_blocks_allChr.bed


## Overlapping GWAS variants with GRBs

After LD block calculation, the user should edit the script run\_targPred.sh by introducing the path to the directory containing the Rscript executable, as well as the path to the targPred code associated with the COREDIR variable. With these options correctly set up, the user can submit the main targPred script by executing the following command:

```
qsub run_targPred.sge
```

In the first step, the script generate grb granges with target gene information in metadata by executing the script `adding_target_genes_to_grbs.R`. 
The  R script `gwas_snps_in_grb.R` is executed next.. This script processes adds LD blocks to gwas assoication, merges  Reported genes and mapped genes into one (GWAS\_GENES). Also, for a reported SNP there are multiple traits reported so this script separates each trait into a separate row. 

**process_gwas_catalog.R**

Script inputs:
  - `data/inputs/GWAS_2026/gwas-catalog-download-associations-alt-full.tsv`
  
Script outputs:
  - `data/processed/gwas_associations_2026.csv`
  - `data/processed/gwas_snps_2026.txt`

**gwas_snps_in_grb.R**
Log file: `processed/logs/gwas_snps_in_grb.txt`

Script inputs:
  - `data/processed/GRBs_with_targets.gr.rds`
  - `data/inputs/GWAS_2026/gwas-efo-trait-mappings.tsv`
  - `data/inputs/GWAS_2026/gwas-catalog-download-associations-alt-full.tsv`
  

Script outputs:
  - `data/outputs/gwas_associations.csv`
  - `data/outputs/loci_in_GRBs.bed`
  - `data/processed/loci_in_GRBs.df.csv`
  - `data_dir/processed/snps_in_GRBs.gr.rds`

**prepare_for_web**

Script inputs:
   - `data/outputs/gwas_associations.csv`
   - `data/processed/loci_in_GRBs.df.csv`
   - `data/inputs/mart_result.df.csv` (optional, if exists)

Script outputs:
   - `data/outputs/GWAS_SNPs_mart_result.csv` (created if biomaRt query is run or if `mart_result.df.csv` is missing)
   - `data/outputs/web_table.df.csv`
   - `data/outputs/targpred<year>.db` (SQLite database with tables: `grb_table`, `loci_table`)



Finally, the script outputs byGRB\_list structure containing information about all genes contained in each GRB and their interactions with enhancers. This structure is required for subsequent plot geeration.

The last step in targPred workflow is plot production. Two types of plots are generated: Gviz plots depicting genomic regions of interest and beanplots showng transcript expression depending on corresponding enhancer activity and empirical p value association.

