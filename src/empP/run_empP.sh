#!/bin/bash
#SBATCH --job-name=empP_parallel
#SBATCH --output=empP_parallel_%j.out
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10

export CH_IMAGE_STORAGE=$HOME/charliecloud

for chr in chr{1..22} chrX; do

    ch-run -b /home/kprskalo/targPred:/workspace --cd /workspace targpred_v6 -- \
        env HOME=/workspace Rscript src/empP/two_tailed_empP_parallel.R $chr
done


# Directory where chromosome outputs are stored
output_dir="/home/kprskalo/targPred/data/processed/empP"


# Header line
header="enh_pos gene transcript_position EP_dist empP_median_oneTailed empP_median_twoTailed E_active eff_size"

for chr in chr{1..22} chrX; do
    (
        echo "Processing $chr"
        chr_dir="${output_dir}/${chr}"
        cat_file="${output_dir}/empP_CAGEenhPermissive_hg19Ens_6.4Mb_${chr}.txt"
        # Concatenate second lines from all files in the chromosome directory, dropping last column
        find "$chr_dir" -maxdepth 1 -type f -exec sed -n '2p' {} \; | awk '{NF--; print}' OFS=' ' > "$cat_file"
        # Add header to the top of the file
        sed -i "1i $header" "$cat_file"
    ) &
done

wait

echo "All chromosome jobs completed."
