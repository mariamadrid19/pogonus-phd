#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=call_snps
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o call_snps.%j.out
#SBATCH --array=1-11

# Array index (0-based)
ID=$((SLURM_ARRAY_TASK_ID - 1))

# Load tools
module load BCFtools/1.12-GCC-10.3.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate variant_tools

# Reference and directories
REFNAME="Pchal_Bar_SW"
REF="/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/Pchalceus_SW.sorted.fasta"
FAI="${REF}.fai"

VCF_DIR="vcfs"
OUT_DIR="final-vcfs"
mkdir -p "$VCF_DIR" "$OUT_DIR"

# Sample names
samples=()
for i in $(seq -w 1 192); do
  samples+=("Pc25Np${i}")
done

# Chromosomes
CHR=(CHR1 CHR2 CHR3 CHR4 CHR5 CHR6 CHR7 CHR8 CHR9 CHR10 CHR11)
names=(1 2 3 4 5 6 7 8 9 10 11)

chr="${CHR[$ID]}"
name="${names[$ID]}"

# Make list of BAMs
ALL_LIST=$(printf " bams/%s.${REFNAME}.filtered.sorted.bam" "${samples[@]}")

echo "====================================="
echo "Chromosome: $chr ($name)"
echo "Output prefix: $REFNAME.chr_$name"
echo "Calling variants..."
echo "====================================="

# Step 1: Call variants
bcftools mpileup -Oz --threads 36 --fasta-ref "$REF" --regions "$chr" $ALL_LIST --annotate FORMAT/DP | bcftools call -m -Oz -f GQ -o "$VCF_DIR/$REFNAME.chr_${name}.vcf.gz"

bcftools index --csi "$VCF_DIR/$REFNAME.chr_${name}.vcf.gz"

# Step 2: Filter
vcftools --gzvcf "$VCF_DIR/$REFNAME.chr_${name}.vcf.gz" \
  --max-missing 0.8 \
  --minQ 30 \
  --maf 0.05 \
  --remove-indels \
  --recode \
  --stdout \
| bgzip > "$VCF_DIR/$REFNAME.chr_${name}.filtered.vcf.gz"

bcftools index --csi "$VCF_DIR/$REFNAME.chr_${name}.filtered.vcf.gz"

# Step 3: Thin
vcftools --gzvcf "$VCF_DIR/$REFNAME.chr_${name}.filtered.vcf.gz" \
  --thin 5000 \
  --recode \
  --stdout \
| bgzip > "$VCF_DIR/$REFNAME.chr_${name}.filtered.thinned.vcf.gz"

bcftools index --csi "$VCF_DIR/$REFNAME.chr_${name}.filtered.thinned.vcf.gz"

# Step 4: Normalize (split multiallelics)
bcftools norm -m -any \
  -Oz \
  -o "$OUT_DIR/$REFNAME.chr_${name}.filtered.thinned.multiSplit.vcf.gz" \
  "$VCF_DIR/$REFNAME.chr_${name}.filtered.thinned.vcf.gz"

bcftools index --csi "$OUT_DIR/$REFNAME.chr_${name}.filtered.thinned.multiSplit.vcf.gz"

echo "Finished chromosome $chr → final output: $OUT_DIR/$REFNAME.chr_${name}.filtered.thinned.multiSplit.vcf.gz"
