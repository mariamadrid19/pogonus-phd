#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --partition=batch_sapphirerapids
#SBATCH --job-name=winpca
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o pca_windows.%A_%a.out
#SBATCH --array=1-11

cd /scratch/leuven/357/vsc35707/winpca/

# env
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate winpca
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

# chromosomes & lengths
chroms=(CHR1 CHR2 CHR3 CHR4 CHR5 CHR6 CHR7 CHR8 CHR9 CHR10 CHR11)
lens=(53807401 53230543 46351446 44639843 65146079 36339123 47275045 35503798 27201578 22848598 49412164)

# 0-based index for arrays
i=$((SLURM_ARRAY_TASK_ID - 1))
chr="${chroms[$i]}"
len="${lens[$i]}"

# VCF number matches array task id (1..11)
num="${SLURM_ARRAY_TASK_ID}"

# cleaned VCF input
clean_vcf="clean-vcfs/${chr}.clean2.vcf.gz"

region="${chr}:1-${len}"
threads="${SLURM_CPUS_PER_TASK:-24}"

echo "========================================================================"
echo "==> Working on ${chr}"
echo "Input VCF:   ${clean_vcf}"
echo "Region:      ${region}"
echo "Threads:     ${threads}"
echo

echo "[*] Running windowed PCA (GT, 200 kb windows, 100 kb step, no MAF filter, ignore PASS) ..."
# Windowed PCA (bp windows as you originally used)
winpca pca "${chr}" "${clean_vcf}" "${region}" -v GT -w 200000 -i 100000 -m 0 --np -t "${threads}"

#echo "[*] Creating chromplot ..."
# chromplot (uses the PREFIX = $chr outputs from the previous step)
winpca chromplot "${chr}" "${region}" -m metadata.fixed.tsv -g habitat -i 5

echo "==> Done: ${chr}"
