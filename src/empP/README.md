# empP Workflow

This folder contains the prerequisite workflow for calculating empirical p-values for enhancer-promoter pairs. Enhancer-gene associations are quantified using FANTOM5 CAGE expression data, with a modified empirical p-value method that tests both activation and repression. For each FANTOM5-detected enhancer, candidate promoter genes within +/-6.4 Mb are evaluated by comparing promoter expression when the enhancer is active (`E+`, CAGE signal > 0) versus inactive (`E0`, CAGE signal = 0), using 10,000 permutations of biosample group labels to estimate two-tailed empirical p-values.

The resulting chromosome-level empP files are used by the core targPred pipeline when building GRB-level summaries and plotting enhancer-promoter relationships.

## Inputs

| Name | Type | Description | Location |
|---|---|---|---|
| `fantom5.enhancers.tpm.pooled.gr.rds` | RDS / GRanges | FANTOM5 enhancer expression object pooled across cell or tissue types. | `data/inputs/` |
| `promoter_cage_tpm_df_fixed_pooled.rds` | RDS / data frame | Promoter CAGE expression matrix used to compare promoter expression across samples. | `data/inputs/` |
| `TSSs.perTranscript.gr.rds` | RDS / GRanges | Transcript-level transcription start site ranges used to define promoter positions. | `data/inputs/` |

## Outputs

| Name | Type | Description | Location |
|---|---|---|---|
| `chr*/<enhancer>_<gene>.txt` | Text files | Per enhancer-promoter pair empirical p-value results produced separately for each chromosome. | `data/processed/empP/` |
| `empP_CAGEenhPermissive_hg19Ens_6.4Mb_chr*.txt` | Text files | Per-chromosome concatenated empP files used by the core targPred pipeline. | `data/processed/empP/` |
| `empP_log.txt` | Log file | Log of the number of enhancer-promoter pairs tested per chromosome. | `data/processed/logs/` |
| `empP_parallel_<jobid>.out` | Scheduler log | Standard output from the empP job script. | Project working directory |

## How to Run

The main entry point is the job script:

```bash
bash src/empP/run_empP.sh
```

The script loops over chromosomes `chr1` through `chr22` and `chrX`, running the chromosome-level R script inside the container:

```bash
Rscript src/empP/two_tailed_empP_parallel.R chr1
```

After chromosome-level pair files are produced, `run_empP.sh` concatenates the per-pair outputs into one empP file per chromosome.

## Main Scripts

- `run_empP.sh`: job script that runs the chromosome-level empP calculation and concatenates per-pair outputs into chromosome-level files.
- `two_tailed_empP_parallel.R`: calculates empirical p-values for enhancer-promoter pairs on one chromosome.
- `empP_utils.R`: helper functions for creating enhancer-promoter combinations, calculating empirical p-values, computing effect sizes, and working with genomic ranges.

## Dependencies / Environment

TODO: add Docker or container image details.

The workflow is currently written to run in a container through `ch-run` and includes scheduler directives for job submission.

Key R packages used by this workflow include:

- `GenomicRanges`
- `dplyr`
- `stringr`
- `stats`
- `parallel`

## Notes

The empP workflow must be completed before the core targPred pipeline steps that require enhancer-promoter empirical p-values. The core workflow expects the chromosome-level empP files to be available under `data/processed/empP/`.
