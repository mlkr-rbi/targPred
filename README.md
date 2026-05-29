# targPred

Genome-wide association studies (GWAS) are a key tool for discovering associations between single nucleotide polymorphisms (SNPs) and phenotypic traits, but identifying the gene affected by a disease-associated SNP remains difficult, especially when the causal gene is distal in genomic distance. `targPred` is an analysis project for linking GWAS loci to candidate target genes through genomic regulatory blocks (GRBs), linkage disequilibrium blocks, and enhancer-promoter evidence.

The project uses the GRB model to infer connections between SNPs or loci and target genes located in the same regulatory block. This approach is motivated by the fact that many non-coding disease-associated variants may act through long-range regulation rather than the nearest gene. The analysis highlights disease traits with a propensity for long-range regulation, including cancer and psychiatric disease, and includes a childhood-obesity-associated locus connected to the distal `BDNF` gene.

The workflow is organized into three main parts: `src/empP` calculates prerequisite enhancer-promoter empirical p-values, `src/core` runs the main targPred pipeline, and `src/analysis` contains downstream analyses based on intermediate or final core outputs.

Main results are available through the targPred website:

```text
http://targpred.irb.hr/
```

Paper: TODO add manuscript/preprint/publication link.

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

The documented workflow focuses on `src/empP`, `src/core`, and the active analysis folders under `src/analysis`.

Large data files, generated outputs, local logs, keys, and local environment folders should not be committed to GitHub unless explicitly intended.
