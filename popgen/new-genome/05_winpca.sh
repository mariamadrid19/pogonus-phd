#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=winpca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o pca_windows.%A_%a.out
#SBATCH --array=1-11

# === Setup ===
cd /scratch/leuven/357/vsc35707/winpca/
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate winpca

# === Define chromosomes and lengths ===
CHRS=(CHR1 CHR2 CHR3 CHR4 CHR5 CHR6 CHR7 CHR8 CHR9 CHR10 CHR11)
CHR_LENGTHS=(53807401 53230543 46351446 44639843 65146079 36339123 47275045 35503798 27201578 22848598 49412164)

# Convert array index (1–11) to 0-based
ID=$((SLURM_ARRAY_TASK_ID - 1))
CHR=${CHRS[$ID]}
LEN=${CHR_LENGTHS[$ID]}

# === Paths ===
# The VCF needs to be split for multiallelic SNPs (could also be imputed with BEAGLE)
VCF="vcfs/Pchal_Bar_SW.chr_${SLURM_ARRAY_TASK_ID}.split.variants.filtered.vcf.gz"

# === Run PCA ===
echo ">>> Running winpca on ${CHR} (${LEN} bp)"
winpca pca ${CHR} "${VCF}" "${CHR}:1-${LEN}" -t 12

echo ">>> Done with ${CHR}"
