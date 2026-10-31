# targPred Core Workflow

This workflow includes calculating LD blocks from GWAS data, assigning GWAS and GRB target genes to SNPs, creating a by-GRB data structure, and creating plots.

The scripts with scheduler directives are to be used in a job scheduling HPC cluster setting. The current job scripts use SLURM-style directives and containerized execution with `ch-run`; older notes may refer to SGE scripts, which can be adapted for usage on other job scheduling systems such as Slurm or PBS.

## Inputs

| Name | Type | Description | Location |
|---|---|---|---|
| `paths.txt` | Text config | Path configuration file used by the core R scripts. Lines beginning with `$` define keywords; the following line contains the value. | `src/core/paths.txt` |
| `gwas-catalog-download-associations-alt-full.tsv` | TSV | GWAS catalog associations file used to parse SNP, trait, and gene information. For the 2026 run, this is downloaded from `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/02/17/`. | `data/inputs/GWAS_<year>/` |
| `gwas-catalog-ancestry.tsv` | TSV | GWAS catalog ancestry file used to retain studies from European ancestry. For the 2026 run, this is downloaded from `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/02/17/`. | `data/inputs/GWAS_<year>/` |
| `gwas-efo-trait-mappings.tsv` | TSV | GWAS-to-EFO trait mappings used when separating multi-trait GWAS rows. For the 2026 run, this is downloaded from `https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/02/17/`. | `data/inputs/GWAS_<year>/` |
| `hg19.canFam3.98.50.processed.threshold.grbs.bed` | BED | Main GRB interval dataset. | `data/inputs/` |
| `hg19.canFam3.tsv` | TSV | GRB target gene annotations. | `data/inputs/` |
| `genes.hg19.ensembl.rds` | RDS | Ensembl gene coordinates used for gene classification in GRBs. | `data/inputs/` |
| `TSSs.perTranscript.gr.rds` | RDS / GRanges | Transcript-level TSS ranges. | `data/inputs/` |
| `fantom5.enhancers.tpm.pooled.gr.rds` | RDS / GRanges | FANTOM5 enhancer expression object. | `data/inputs/` |
| `promoter_cage_tpm_df_fixed_pooled.rds` | RDS / data frame | Promoter CAGE expression matrix. | `data/inputs/` |
| `empP_CAGEenhPermissive_hg19Ens_6.4Mb_chr*.txt` | Text files | Chromosome-level empP files produced by the prerequisite empP workflow. | `data/processed/<empP_dir>/` |
| PLINK / 1000 Genomes files | PLINK files | Reference genotype data used to calculate LD blocks. `run_PLINK.sh` can download and prepare these files. | `data/plink_dir/` |
| `gwas_coordinates.tsv` | TSV | Local 1000 Genomes Phase 3 GRCh37 SNP coordinates covering every LD lead ID. | `data/inputs/1000genomes_phase3_grch37/` |

## Outputs

| Name | Type | Description | Location |
|---|---|---|---|
| `gwas_associations_<year>.csv` | CSV | Filtered GWAS associations after ancestry and SNP filtering. | `data/processed/` |
| `gwas_snps_<year>.txt` | Text file | Unique GWAS SNP list used for LD block calculation. | `data/processed/` |
| `EUR_LD_blocks_allChr.bed` | BED | Combined LD block file for all chromosomes. | `data/processed/` |
| `GRBs_with_targets.gr.rds` | RDS / GRanges | GRB intervals with target gene metadata. | `data/processed/` |
| `GRBs_with_targets.tsv` | TSV | Table version of GRBs with target gene metadata. | `data/processed/` |
| `gwas_associations.csv` | CSV | GWAS associations joined to LD block information and normalized trait/gene columns. | `data/outputs/` |
| `loci_in_GRBs.bed` | BED | GWAS loci overlapping GRBs for genome-browser visualization. | `data/outputs/` |
| `loci_in_GRBs.df.csv` | CSV | Data frame of loci in GRBs for downstream analyses. | `data/processed/` |
| `snps_in_GRBs.gr.rds` | RDS / GRanges | GRanges object containing SNPs that fall in GRBs. | `data/processed/` |
| `web_table.csv` | CSV | Main table prepared for downstream/web use. | `data/outputs/` |
| `targpred<year>.db` | SQLite database | SQLite database with `grb_table` and `loci_table`. | `data/outputs/` |
| `<byGRB_name>` | RDS / list | by-GRB list containing genes in each GRB and enhancer-promoter information. | `data/processed/` |
| `gene_model_hg19.gr.rds` | RDS / GRanges | Transcript gene model from GRCh37 Ensembl. | `data/processed/` |
| Beanplots | PNG | Transcript expression plots by enhancer activity and empirical p-value association. | `data/plots/expression_pairs/` |
| Gviz plots | PNG | Genomic-region plots for loci and GRBs. | `data/plots/gviz_plots/` |
| Core logs | Text files | Processing logs for GWAS, PLINK, GRB overlap, and by-GRB creation. | `data/processed/logs/` |

## Path Configuration

Before running the core workflow, edit `src/core/paths.txt` with the required paths and filenames for targPred to function properly. The file contains keywords and values. The lines containing keywords start with the `$` sign, while the following line must contain the appropriate value for that keyword. Any line beginning with a `#` sign will be interpreted as a comment and ignored.

Meaning of the keywords:

- `$core`: path to your targPred code directory where the scripts are located.
- `$data`: path to the directory containing data, organized as follows:
  1. `inputs`: directory containing `EFO_terms.rds`, GWAS catalog files, `hg19.canFam3.98.50.processed.threshold.grbs.bed`, and `hg19.canFam3.tsv`.
  2. `processed`: directory containing empP files, `genes.hg19.ensembl.rds`, `TSSs.perTranscript.gr.rds`, `promoter_cage_tpm_df_fixed_pooled.rds`, and `fantom5.enhancers.tpm.pooled.gr.rds`. Some intermediary files are also stored here.
  3. `outputs`: resulting files such as the web table and SQLite database are stored here.
- `$work`: path to your working directory.
- `$empP_dir`: directory name with your empP files.
- `$empP`: file prefix of empP files.
- `$year`: year of the GWAS catalog edition.
- `$LDblock_path`: full path leading to the bed file with LD blocks, resulting from the PLINK analysis. The expected file name is `EUR_LD_blocks_allChr.bed`.
- `$byGRB_name`: filename used to save the by-GRB list.

## Dependencies / Environment

TODO: add Docker or container image details.

External tools used by the core workflow include:

- PLINK and PLINK2
- bedtools
- wget and unzip
- R / Rscript
- a job scheduler environment for selected shell scripts

## Downloading GWAS Data and Calculating LD Blocks Using 1000 Genomes Data

To process the GWAS catalog, use the `process_gwas_catalog.R` script. This script parses the GWAS catalog release, relevant SNP and trait information, and formats the data for downstream analysis.

For the 2026 run, `run_PLINK.sh` downloads the GWAS release files from:

```text
https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/02/17/
```

The files from this release are stored under `data/inputs/GWAS_2026/` and include the associations file, ancestry file, and EFO trait mapping file.

The script `run_PLINK.sh` follows the workflow outlined in `https://cran.r-project.org/web/packages/snpsettest/vignettes/reference_1000Genomes.html`. It downloads the hg19 data from `https://www.cog-genomics.org/plink/2.0/resources#1kg_phase3`, selects the desired population, transforms the data into PLINK1 format, and calculates the linkage disequilibrium blocks.

The user should add the absolute or relative paths to directories containing PLINK1 and PLINK2 executables, as well as the path to the text file containing the list of SNPs for the LD block calculation. The user can change the PLINK settings for LD block calculations directly in the script.

The script can afterwards be submitted to the queueing system, or run directly if the same environment is available:

```bash
bash src/core/run_PLINK.sh
```

The bed file with the linkage disequilibrium blocks can be found in the newly created directory named `ld_blocks`. The script produces a file with LD blocks for each chromosome separately, as well as their concatenated version. The bed file with all combined chromosomes is stored in `data/processed/EUR_LD_blocks_allChr.bed`.

## Overlapping GWAS Variants With GRBs

After LD block calculation, the user should edit `run_targPred.sh` by introducing the path to the directory containing the `Rscript` executable, as well as the path to the targPred code associated with the `COREDIR` variable. With these options correctly set up, the user can submit the main targPred script to the queueing system, or run it directly if the same environment is available:

```bash
bash src/core/run_targPred.sh
```

In the first step, the script generates GRB GRanges with target gene information in metadata by executing `adding_target_genes_to_grbs.R`.

The R script `gwas_snps_in_grb.R` is executed next. This script adds LD blocks to GWAS associations and merges reported genes and mapped genes into one `GWAS_GENES` column. Also, for a reported SNP there are multiple traits reported, so this script separates each trait into a separate row.

Finally, the workflow outputs a by-GRB list structure containing information about all genes contained in each GRB and their interactions with enhancers. This structure is required for subsequent plot generation.

The last step in the targPred workflow is plot production. Two types of plots are generated: Gviz plots depicting genomic regions of interest, and beanplots showing transcript expression depending on corresponding enhancer activity and empirical p-value association.

## Run Order and Script Outputs

### `process_gwas_catalog.R`

Script inputs:

- `data/inputs/GWAS_<year>/gwas-catalog-download-associations-alt-full.tsv`
- `data/inputs/GWAS_<year>/gwas-catalog-ancestry.tsv`

Script outputs:

- `data/processed/gwas_associations_<year>.csv`
- `data/processed/gwas_snps_<year>.txt`
- `data/processed/logs/process_gwas_catalog.txt`

### `run_PLINK.sh`

Script inputs:

- `data/processed/gwas_snps_<year>.txt`
- PLINK / PLINK2 executables
- 1000 Genomes reference files downloaded or stored under `data/plink_dir/`

Script outputs:

- per-chromosome LD files under `data/plink_dir/ld_blocks/`
- `data/processed/EUR_LD_blocks_allChr.bed`
- `data/processed/logs/plink_processing.log`

### `adding_target_genes_to_grbs.R`

Script inputs:

- `data/inputs/hg19.canFam3.98.50.processed.threshold.grbs.bed`
- `data/inputs/hg19.canFam3.tsv`

Script outputs:

- `data/processed/GRBs_with_targets.gr.rds`
- `data/processed/GRBs_with_targets.tsv`
- `data/processed/logs/targPred_run_log.txt`

### `gwas_snps_in_grb.R`

Log file: `data/processed/logs/gwas_snps_in_grb.txt`

SNP positions come from the local 1000 Genomes Phase 3 GRCh37 lookup at `data/inputs/1000genomes_phase3_grch37/gwas_coordinates.tsv` (columns `refsnp_id`, `chr_name`, `chrom_start`; 1-based positions). Include every lead ID from the LD-block BED when preparing this lookup, even IDs absent from the current GWAS SNP list. Before writing results, the script checks that each LD lead has one valid position on the expected primary chromosome (1–22 or X). Missing IDs, conflicting positions, or chromosome disagreements stop processing. The web table takes its chromosome from the LD locus and joins positions by both rsID and chromosome. No BioMart query or liftover is required.

Script inputs:

- `data/processed/GRBs_with_targets.gr.rds`
- `data/processed/gwas_associations_<year>.csv`
- `data/inputs/GWAS_<year>/gwas-efo-trait-mappings.tsv`
- `data/processed/EUR_LD_blocks_allChr.bed`
- `data/inputs/1000genomes_phase3_grch37/gwas_coordinates.tsv`

Script outputs:

- `data/outputs/gwas_associations.csv`
- `data/outputs/loci_in_GRBs.bed`
- `data/processed/loci_in_GRBs.df.csv`
- `data/processed/snps_in_GRBs.gr.rds`
- `data/outputs/web_table.csv`
- `data/outputs/targpred<year>.db`

### `create_byGRB_list.R`

Script inputs:

- `data/processed/GRBs_with_targets.gr.rds`
- `data/processed/loci_in_GRBs.df.csv`
- chromosome-level empP files under `data/processed/<empP_dir>/`
- `data/inputs/genes.hg19.ensembl.rds`
- `data/inputs/TSSs.perTranscript.gr.rds`
- `data/inputs/fantom5.enhancers.tpm.pooled.gr.rds`
- `data/inputs/promoter_cage_tpm_df_fixed_pooled.rds`

Script outputs:

- `data/processed/<byGRB_name>`
- `data/processed/gene_model_hg19.gr.rds`
- `data/processed/logs/byGRB_list_processing.txt`

### `locus_summary_plots.R`

Script inputs:

- `data/processed/<byGRB_name>`
- `data/processed/loci_in_GRBs.df.csv`
- expression and GRB input files listed above

Script outputs:

- beanplots under `data/plots/expression_pairs/`

### `locus_summary_plots_gviz.R`

Script inputs:

- `data/processed/<byGRB_name>`
- `data/processed/loci_in_GRBs.df.csv`
- `data/processed/snps_in_GRBs.gr.rds`
- expression, GRB, HCNE, and TAD-track input files listed above

Script outputs:

- Gviz plots under `data/plots/gviz_plots/`
- `data/plots/gviz_plots/gviz_status.tsv`: one row per trait/GRB panel with status `written`, `skipped` (the plotting helper found no usable loci or enhancers) or `error` and the error message
- `data/cache/gviz/`: GRCh37 Mart connection, gene tracks per GRB region, ideograms and hg19 cytobands (`hg19_cytoBandIdeo.txt.gz`, downloaded from UCSC on first use)

The script runs traits in parallel with `SLURM_CPUS_PER_TASK` cores (one if unset) and stops with an error after writing the status file if any panel failed. Network access and caching are handled by `gviz_cache.R`:

- Gene tracks come from the GRCh37 Ensembl archive (`grch37.ensembl.org`). The connection is only opened on a cache miss, so a run with a populated `data/cache/gviz/` does not contact Ensembl or UCSC for track construction. Keep the cache with the data; do not reuse it for a different assembly or gene-annotation release.
- Gviz queries Ensembl again while drawing whenever a panel extends beyond the cached gene track, which is the case for most panels. Panels failing with transient Ensembl errors (HTTP 5xx/405, timeouts, connection failures) are retried up to five times with waits of 1, 2, 4 and 8 minutes. Other errors are not retried.
- Each worker uses its own biomaRt query cache under the local temporary directory. The default `~/.cache/biomaRt` serialises parallel workers on its sqlite lock when the home directory is on a network file system.
- biomaRt's hard-coded 10-second request timeout is raised to 180 seconds.

Package versions: Gviz 1.42.1 requires knitr >= 1.43 at load time, and the `collect(Inf)` call in BiocFileCache 2.6.1 (used by biomaRt 2.54.1) fails with dbplyr 2.5.0 but works with dbplyr 2.3.4. The shared `targPred` conda environment has knitr 1.40 and dbplyr 2.5.0, so install knitr 1.46 and dbplyr 2.3.4 into a separate library and point `R_LIBS_USER` to it (as in `runs/coordinate_refresh_20261027/R-library`).
