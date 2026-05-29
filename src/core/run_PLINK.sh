#!/bin/bash
#SBATCH --job-name=targPred_job
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

export CH_IMAGE_STORAGE=$HOME/charliecloud

GWAS_DIR=data/inputs/GWAS_2026
GWAS_RELEASE=2026/02/17/


SNPLIST=data/processed/gwas_snps_2026.txt

# Create directory if it doesn't exist and it will store all the PLINK files
PLINKDIR=data/plink_dir
mkdir -p $PLINKDIR


# Define where plink and plink2 is located and add it to PATH
export PLINK="data/plink_dir/plink"
export PATH=$PATH:$PLINK
    
export PLINK2="data/plink_dir/plink2"
export PATH=$PATH:$(dirname $PLINK2)

POP="EUR"

# Create a log file in data/processed/logs
LOG_DIR="data/processed/logs"
mkdir -p $LOG_DIR
LOG_FILE="$LOG_DIR/plink_processing.log"

# Define the output file for the final BED file with LD blocks
OUTPUT_BED="data/processed/${POP}_LD_blocks_allChr.bed"

###### Download GWAS catalog data #######
# this is showcased on GWAS 2026 https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/02/17/

mkdir -p $GWAS_DIR

# Downloading associations with ontology annotation
wget -O $GWAS_DIR/gwas-catalog-associations_ontology-annotated-full.zip "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/$GWAS_RELEASE/gwas-catalog-associations_ontology-annotated-full.zip"

unzip $GWAS_DIR/gwas-catalog-associations_ontology-annotated-full.zip -d $GWAS_DIR

# Check number of lines in gwas-catalog-associations_ontology-annotated-full.tsv
NUM_ASSOCIATIONS=$(wc -l < $GWAS_DIR/gwas-catalog-associations_ontology-annotated-full.tsv)
echo "Number of associations in the GWAS catalog: $NUM_ASSOCIATIONS" >> $LOG_FILE

# Downloading ancestry information
wget -O $GWAS_DIR/gwas-catalog-ancestry.tsv "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/$GWAS_RELEASE/gwas-catalog-ancestry.tsv"

# Downloading EFO traits
wget -O $GWAS_DIR/gwas-efo-trait-mappings.tsv "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/$GWAS_RELEASE/gwas-efo-trait-mappings.tsv"

# Clean up GWAS data to select only European ancestry and prepare SNP list for LD calculation
ch-run -b /home/kprskalo/targPred:/workspace --cd /workspace targpred_v6 -- env HOME=/workspace Rscript src/core/process_gwas_catalog.R $GWAS_DIR

# Add number of SNPs to log file
NUM_SNPS=$(wc -l < $SNPLIST)
echo "Number of SNPs in the list: $NUM_SNPS" >> $LOG_FILE

echo "-----starting to download dataset----"

# The links in here may be changed in future
# "-O" to specify output file name
wget -O $PLINKDIR/all_phase3.psam "https://www.dropbox.com/s/6ppo144ikdzery5/phase3_corrected.psam?dl=1"
wget -O $PLINKDIR/all_phase3.pgen.zst "https://www.dropbox.com/s/y6ytfoybz48dc0u/all_phase3.pgen.zst?dl=1"
wget -O $PLINKDIR/all_phase3.pvar.zst "https://www.dropbox.com/s/odlexvo8fummcvt/all_phase3.pvar.zst?dl=1"

# Decompress pgen.zst to pgen
$PLINK2 --zst-decompress $PLINKDIR/all_phase3.pgen.zst > $PLINKDIR/all_phase3.pgen

echo "-----finished downloading the dataset----"

echo "-----extracting the desired population----"


$PLINK2 --pfile $PLINKDIR/all_phase3 vzs --chr 1-22,X --output-chr 26 --max-alleles 2 --rm-dup exclude-mismatch --set-missing-var-ids '@_#_$1_$2' --make-pgen --allow-extra-chr --new-id-max-allele-len 500 --out $PLINKDIR/all_phase3_autosomes_X --memory 200000

# Prepare sub-population filter file
# here $6 == CEU or IBS or whatever
# but $5 == EUR
awk 'NR == 1 || $5 == "EUR" {print $1}' $PLINKDIR/all_phase3_autosomes_X.psam > $PLINKDIR/${POP}_1kg_samples_autosomes_X.txt


# Generate sub-population fileset for all
$PLINK2 --pfile $PLINKDIR/all_phase3_autosomes_X --keep $PLINKDIR/${POP}_1kg_samples_autosomes_X.txt --make-pgen --allow-extra-chr --new-id-max-allele-len 500 --out $PLINKDIR/${POP}_phase3_autosomes_X 

echo "-----finished extracting the desired population----"

echo "-----separating the data by chromosome and running LD block calculations----"


# Prepare plink2 output for plink (all chromosomes)
$PLINK2 --pfile $PLINKDIR/${POP}_phase3_autosomes_X --maf 0.005 --make-bed --allow-extra-chr --new-id-max-allele-len 500 --out $PLINKDIR/${POP}_phase3_autosomes_X

mkdir -p $PLINKDIR/plink_files $PLINKDIR/ld_blocks

for i in {1..22} X
do 
    NAME=${POP}_phase3_chr${i}
    $PLINK2 --bfile $PLINKDIR/${POP}_phase3_autosomes_X --chr $i --make-bed --allow-extra-chr --new-id-max-allele-len 500 --out $PLINKDIR/plink_files/${NAME}
    $PLINK --bfile $PLINKDIR/plink_files/${NAME} --r2 --ld-snp-list ${SNPLIST} --ld-window-kb 500 --ld-window 99999 --ld-window-r2 0.6 --noweb --allow-extra-chr --out $PLINKDIR/ld_blocks/${NAME}
done

echo "-----LD block calculation finished----"

echo "-------creating BED files with LD blocks for targPred--------"

# Run R script 
ch-run -b /home/kprskalo/targPred:/workspace --cd /workspace targpred_v6 -- env HOME=/workspace Rscript src/core/make_bed.R $POP $PLINKDIR $OUTPUT_BED


# Add number of LD blocks to log file
NUM_BLOCKS=$(wc -l < data/processed/${POP}_LD_blocks_allChr.bed)
echo "Number of LD blocks identified: $NUM_BLOCKS" >> $LOG_FILE

# Check how many of the GWAS SNPs are not in the LD blocks and add that to log file
echo "GWAS SNPs not in LD blocks:" >> $LOG_FILE
grep -Fvxf <(cut -f4 data/processed/${POP}_LD_blocks_allChr.bed) $SNPLIST >> $LOG_FILE


echo "---------job finished-------------"
