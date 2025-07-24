#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=thin_snps
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=16G
#SBATCH --output=thin_snps_%A_%a.out
#SBATCH --error=thin_snps_%A_%a.err
#SBATCH --array=0-4
#SBATCH -A lp_svbelleghem

# Chromosome mapping: index to actual chromosome numbers
chromosomes=(2 4 5 6 8)
chr=${chromosomes[$SLURM_ARRAY_TASK_ID]}

# Directories
INPUT_DIR="filtered-no-missing"
OUTPUT_DIR="filtered-thinned"
mkdir -p "$OUTPUT_DIR"

input_vcf="${INPUT_DIR}/Pchal_Bar_SW.chr_${chr}.filtered.nomissing.vcf.gz"
output_vcf="${OUTPUT_DIR}/Pchal_Bar_SW.chr_${chr}.filtered.nomissing.thin500.vcf.gz"

echo "Processing CHR${chr}..."

# Check for .csi index
if [ ! -f "${input_vcf}.csi" ]; then
  echo "Index (.csi) not found for $input_vcf, creating it..."
  bcftools index "$input_vcf"
fi

# Perform SNP thinning (min distance 500bp)
bcftools view -v snps "$input_vcf" | \
  bcftools +prune -Ou --distance 500 | \
  bcftools view -Oz -o "$output_vcf"

bcftools index "$output_vcf"

echo "Done thinning CHR${chr}. Output saved to $output_vcf"
