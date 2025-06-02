#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name spain_snps
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o spain_snps.%j.out
#SBATCH --array=1-263

# Index starts at 0
ID=$((SLURM_ARRAY_TASK_ID - 1))

# Load required modules
module load SAMtools/0.1.20-GCC-12.3.0
module load Python/3.12.3-GCCcore-13.3.0
module load tabixpp/1.1.0-GCC-10.3.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

# Define reference
REFNAME=P_chalceus_REF1
REF=/scratch/leuven/357/vsc35707/popgen/${REFNAME}.fa
FAI=${REF}.fai

# Read scaffold names from the .fai file into an array
readarray -t scaffolds < <(cut -f1 $FAI)

# Get current scaffold
CHR=${scaffolds[$ID]}

# Sample IDs
samples=(GC136107 GC136108 GC136109 GC136110 GC136111 GC136112 GC136113 GC136114 GC136115 GC136116)

# Generate list of BAM files
ALL_LIST=""
for SAMPLE in "${samples[@]}"; do
    ALL_LIST+=" $SAMPLE.${REFNAME}.filtered.sorted.nd.bam"
done

# Output file stem
OUTBASE="Pogonus_${REFNAME}.${CHR}"

# Check validity
if [[ -z "$CHR" ]]; then
    echo "Error: scaffold ID is empty!"
    exit 1
fi

# Variant calling
bcftools mpileup -Oz --threads 20 -f "$REF" $ALL_LIST -r "$CHR" | bcftools call -m -Oz -o ${OUTBASE}.vcf.gz

# Filtering with vcftools
vcftools --gzvcf "${VCF_CALL}.vcf.gz" --recode --remove-indels --minQ 30 --max-missing 0.25 --stdout | bgzip > "${VCF_CALL}.filt.bi.vcf.gz"

# Parsing VCF with filtering rules
/data/leuven/357/vsc35707/miniconda3/bin/python3.11 parseVCF.py \
  --gtf flag=GQ   min=30   gtTypes=Het \
  --gtf flag=GQ   min=30   gtTypes=HomAlt \
  --gtf flag=DP   min=10 \
  --skipIndels \
  -i "${VCF_CALL}.filt.bi.vcf.gz" \
  | gzip > "${VCF_CALL}.calls.filt.bi.vcf.gz"

# Striping BAM suffix from SNP IDs
CALLS_H="Pogonus_${REFNAME}_chr_${CHRNAME}.H.calls"
zcat "${VCF_CALL}.calls.filt.bi.vcf.gz" | sed 's/\.filtered\.sorted\.nd\.bam//g' | bgzip -c > "${CALLS_H}.gz"

echo "Done chr ${CHRNAME}."
