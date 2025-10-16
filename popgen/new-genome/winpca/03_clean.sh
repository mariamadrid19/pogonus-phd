#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --partition=batch_sapphirerapids
#SBATCH --job-name=clean2vcf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=04:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o clean2.%A_%a.out
#SBATCH --array=1-11

# --- config ---
IN_DIR="new-vcfs"
OUT_DIR="clean-vcfs"
DROP_MONO="${DROP_MONO:-1}"  # 1 = drop monomorphic; 0 = keep
THREADS="${SLURM_CPUS_PER_TASK:-4}"

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate winpca
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

# array index -> chromosome number
i="${SLURM_ARRAY_TASK_ID}"
chr="CHR${i}"

in_vcf="${IN_DIR}/${chr}.subset.vcf.gz"
out_vcf="${OUT_DIR}/${chr}.clean2.vcf.gz"

echo "==> ${chr}"
echo " Input : ${in_vcf}"
echo " Output: ${out_vcf}"
echo " Drop monomorphic: ${DROP_MONO}"
echo

# sanity check
if [[ ! -s "${in_vcf}" ]]; then
  echo "ERROR: input VCF not found or empty: ${in_vcf}" >&2
  exit 1
fi

# Choose filter expression
# Always: remove sites where all samples are missing (AN==0 or F_MISSING==1)
# If DROP_MONO=1 also remove monomorphic (AC==0 or AC==AN)
if [[ "${DROP_MONO}" -eq 1 ]]; then
  expr='AN>0 && AC>0 && AC<AN && F_MISSING<1'
else
  expr='AN>0 && F_MISSING<1'
fi

# Clean: keep only biallelic SNPs; annotate AC/AN/F_MISSING; filter; bgzip+index
echo "[*] Cleaning -> drop all-missing $( [[ ${DROP_MONO} -eq 1 ]] && echo 'and monomorphic' || echo '(keep monomorphic)' ) ..."
bcftools view -m2 -M2 -v snps -Ou "${in_vcf}" \
| bcftools +fill-tags -Ou -- -t AC,AN,F_MISSING \
| bcftools view -i "${expr}" -Oz -o "${out_vcf}"

bcftools index -f "${out_vcf}"

# Post-clean checks
echo "[*] Counting sites after..."
after=$(bcftools index -n "${out_vcf}")
echo "    Sites after : ${after}"

echo "[*] Verifying no all-missing records remain..."
left=$(bcftools +fill-tags "${out_vcf}" -- -t AN,F_MISSING \
       | bcftools view -H -i 'AN==0 || F_MISSING==1' | wc -l)
echo "    All-missing remaining: ${left}"

echo "Done ${chr}"
