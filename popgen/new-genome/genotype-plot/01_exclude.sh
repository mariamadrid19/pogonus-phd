#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=exclude
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --output=exclude_%A_%a.out
#SBATCH --error=exclude_%A_%a.err
#SBATCH --array=0-4
#SBATCH -A lp_svbelleghem

# List of raw sample IDs (without full path or suffix)
sample_ids=(GC129399 GC129394 GC129419 GC129420 GC129439 GC129440 GC136105 GC136106)

# Format sample names to match VCF header format
exclude_samples=()
for id in "${sample_ids[@]}"; do
  exclude_samples+=("bams/${id}.Pchal_Bar_SW.filtered.sorted.dedup.bam")
done

# Chromosomes to process
chromosomes=(2 4 5 6 8)

# Directories
INPUT_DIR="."
OUTPUT_DIR="filtered-no-missing"
mkdir -p "$OUTPUT_DIR"

# Exclusion file
EXCLUDE_FILE="samples_to_exclude.txt"
if [ ! -f "$EXCLUDE_FILE" ]; then
  printf "%s\n" "${exclude_samples[@]}" > "$EXCLUDE_FILE"
fi

# Get current chromosome from SLURM array index
chr="${chromosomes[$SLURM_ARRAY_TASK_ID]}"
input_vcf="${INPUT_DIR}/Pchal_Bar_SW.chr_${chr}.filtered.vcf.gz"
vcf_index="${input_vcf}.csi"
output_vcf="${OUTPUT_DIR}/Pchal_Bar_SW.chr_${chr}.filtered.nomissing.vcf.gz"

# Check if input VCF index exists
if [ ! -f "$vcf_index" ]; then
  echo "Index not found for $input_vcf. Indexing now..."
  bcftools index "$input_vcf"
else
  echo "Index found for $input_vcf. Skipping indexing."
fi

# Filter VCF
echo "Filtering chromosome $chr..."

bcftools view \
  -S ^"$EXCLUDE_FILE" \
  -i 'COUNT(GT="mis")==0' \
  -Oz -o "$output_vcf" \
  "$input_vcf"

bcftools index "$output_vcf"

echo "Finished chromosome $chr"
