#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=subset_vcfs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=04:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o subset_vcfs.%j.out

# ============================================================
# Subset each chromosome VCF to include only samples listed in samples.txt
# Creates a new directory 'new-vcfs' with 11 filtered VCFs (CHR1–CHR11)
# ============================================================

cd /scratch/leuven/357/vsc35707/winpca/

# Activate conda if needed
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate winpca

# Input / output directories
VCF_DIR="vcfs"
OUT_DIR="new-vcfs"
SAMPLES="samples.txt"

# Create output directory
mkdir -p "$OUT_DIR"

echo "=========================================================="
echo " Subsetting VCFs to samples in $SAMPLES"
echo " Output directory: $OUT_DIR"
echo "=========================================================="

# Loop over 11 chromosome VCFs
for in_vcf in ${VCF_DIR}/Pchal_Bar_SW.chr_*.split.biallelic.vcf.gz; do
    base=$(basename "$in_vcf")
    out_vcf="${OUT_DIR}/${base/.split.biallelic.vcf.gz/.subset.vcf.gz}"

    echo "→ Processing: $base"
    echo "  Output:     $(basename "$out_vcf")"

    # Keep only listed samples (ignore missing ones)
    bcftools view -S "$SAMPLES" --force-samples -Oz -o "$out_vcf" "$in_vcf"

    # Index the new file
    tabix -f -p vcf "$out_vcf"

    # Report number of samples kept
    count=$(bcftools query -l "$out_vcf" | wc -l)
    echo "  Samples kept: $count"
    echo
done

echo "All done. Filtered files saved in: $OUT_DIR"
