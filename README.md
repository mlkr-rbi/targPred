# targPred

`targPred` is an analysis project for linking GWAS loci to candidate target genes through genomic regulatory blocks (GRBs), linkage disequilibrium blocks, and enhancer-promoter evidence. The project asks whether disease- or trait-associated variants that fall within GRBs can improve target gene prioritization beyond standard GWAS gene annotations.

The workflow is organized into three main parts: `src/empP` calculates prerequisite enhancer-promoter empirical p-values, `src/core` runs the main targPred pipeline, and `src/analysis` contains downstream analyses based on intermediate or final core outputs.

## Inputs and Outputs

Detailed inputs and outputs are documented in the subfolder README files:

- `src/empP/README.md`: prerequisite enhancer-promoter empirical p-value workflow.
- `src/core/README.md`: main targPred workflow, including GWAS processing, LD block generation, GRB overlap, target gene assignment, and plotting.
- `src/analysis/analysis_notes.md`: downstream validation and focused follow-up analyses.

At a high level, required inputs are stored under `data/inputs/`, intermediate files under `data/processed/`, final tables/databases under `data/outputs/`, analysis-specific results under `data/analysis/`, and figures under `data/plots/`.

## Repository Organization

- `src/empP/`: prerequisite code for calculating empirical p-values for enhancer-promoter pairs.
- `src/core/`: main targPred pipeline.
- `src/analysis/`: downstream analyses based on intermediate or final outputs from the core pipeline.
- `src/scripts/`: helper scripts. Currently, `src/scripts/make_bed.R` is called by the core PLINK workflow.
- `src/trash/`: archived or exploratory scripts that are not part of the current documented workflow.
- `data/inputs/`: required input data.
- `data/processed/`: intermediate files produced during the workflow.
- `data/outputs/`: main output tables and database files.
- `data/analysis/`: analysis-specific inputs, outputs, and results.
- `data/plots/`: generated plots.

## Dependencies / Environment

TODO: add Docker or container image details.

The workflow currently uses containerized execution in several job scripts. Some scripts are intended to be submitted to a job scheduler and include scheduler directives.

External tools used by the workflow include:

- PLINK and PLINK2
- bedtools
- wget and unzip
- R / Rscript
- a job scheduler environment for selected shell scripts

## How to Run the Analysis

Run the workflow in this order:

1. Prepare required input files under `data/inputs/`.
2. Run the `src/empP` prerequisite workflow to calculate enhancer-promoter empirical p-values.
3. Run the `src/core` workflow to process GWAS data, calculate LD blocks, overlap loci with GRBs, assign target genes, and generate core outputs.
4. Run selected downstream analyses from `src/analysis`.

See the subfolder documentation for specific commands, expected inputs, and expected outputs.

## Main Workflow Stages

1. `src/empP/`: calculates empirical p-values for enhancer-promoter pairs.
2. `src/core/run_PLINK.sh`: job script for GWAS catalog processing, 1000 Genomes/PLINK setup, LD block calculation, and LD block BED creation.
3. `src/core/run_targPred.sh`: job script for running the main targPred stages inside the container.
4. `src/core/`: R scripts for GRB target gene annotation, GWAS/GRB overlap, by-GRB object creation, and plotting.
5. `src/analysis/`: downstream validation, crude-number summaries, and LRR phenotype analyses.

## Notes

Some scripts are exploratory, archived, or retained for provenance. The documented workflow should focus on `src/empP`, `src/core`, and the active analysis folders under `src/analysis`.

Large data files, generated outputs, local logs, keys, and local environment folders should not be committed to GitHub unless explicitly intended.
