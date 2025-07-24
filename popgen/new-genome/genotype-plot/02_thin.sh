#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=thin_snps
#SBATCH --array=0-4
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=8G
#SBATCH --output=thin_snps_%A_%a.out
#SBATCH --error=thin_snps_%A_%a.err
#SBATCH -A lp_svbelleghem

conda activate variant_tools

chromosomes=(2 4 5 6 8)
chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

INPUT_DIR="filtered-no-missing"
OUTPUT_DIR="filtered-thinned"

input_vcf="${INPUT_DIR}/Pchal_Bar_SW.chr_${chr}.filtered.nomissing.vcf.gz"
output_prefix="${OUTPUT_DIR}/Pchal_Bar_SW.chr_${chr}.filtered.nomissing.thin500"

echo "Thinning CHR${chr}..."
vcftools --thin 500 --gzvcf "$input_vcf" --recode --stdout --remove-indels | bgzip > "${output_prefix}.vcf.gz"

bcftools index --csi "${output_prefix}.vcf.gz"

echo "Finished thinning CHR${chr} → saved to ${output_prefix}.vcf.gz"
