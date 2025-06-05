#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name spain_snps
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o spain_snps.%j.out
#SBATCH --array=1-131

# Index starts at 0
ID=$((SLURM_ARRAY_TASK_ID - 1))

# Load required modules
module load SAMtools/1.16.1-GCC-11.3.0
module load Python/3.12.3-GCCcore-13.3.0
module load tabixpp/1.1.0-GCC-10.3.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

# Reference files
REFNAME="P_chalceus_REF1"
REF="/scratch/leuven/357/vsc35707/popgen/${REFNAME}.fa"
FAI="${REF}.fai"

# Output directories
CALLS_DIR="calls_H"
VCF_DIR="vcfs_H"
mkdir -p "$CALLS_DIR" "$VCF_DIR"

# Get scaffold name from .fai
readarray -t scaffolds < <(cut -f1 "$FAI")
CHRNAME="${scaffolds[$ID]}"

# Sample IDs
samples=(GC136107 GC136108 GC136109 GC136110 GC136111 GC136112 GC136113 GC136114 GC136115 GC136116)

# Build BAM file list from 'bams/' directory
ALL_LIST=""
for SAMPLE in "${samples[@]}"; do
    ALL_LIST+=" bams/${SAMPLE}.${REFNAME}.filtered.sorted.dedup.bam"
done

# Output VCF stem
VCF_CALL="${VCF_DIR}/Pogonus_${REFNAME}.${CHRNAME}"

# Exit early if scaffold name is empty
if [[ -z "$CHRNAME" ]]; then
    echo "Error: empty scaffold name at index $ID"
    exit 1
fi

# Step 1: Call variants with FORMAT annotations
bcftools mpileup -Oz --threads 36 -f "$REF" $ALL_LIST -r "$CHRNAME" -a FORMAT/DP | bcftools call --multiallelic-caller -Oz -f GQ -o "${VCF_CALL}.vcf.gz"

# Index VCF
bcftools index "${VCF_CALL}.vcf.gz"

# Step 2: Filter with vcftools
vcftools --gzvcf "${VCF_CALL}.vcf.gz" --recode --remove-indels --stdout | bgzip > "${VCF_CALL}.filt.bi.vcf.gz"

# Step 3: Parse VCF with custom script
python parseVCF.py \
  --gtf flag=GQ   min=20   gtTypes=Het \
  --gtf flag=GQ   min=20   gtTypes=HomAlt \
  --gtf flag=DP   min=10 \
  --skipIndels \
  -i "${VCF_CALL}.filt.bi.vcf.gz" \
  | gzip > "${VCF_CALL}.calls.gz"

# Step 4: Strip BAM suffix from SNP IDs and save in calls_H/
CALLS_H="${CALLS_DIR}/Pogonus_${REFNAME}_chr_${CHRNAME}.H.calls.gz"
zcat "${VCF_CALL}.calls.gz" \
    | sed 's/\.filtered\.sorted\.dedup\.bam//g' \
    | bgzip -c > "$CALLS_H"

# Step 5: Cleanup only if H.calls file is non-empty
if [[ -s "$CALLS_H" ]]; then
    echo "Successfully created $CALLS_H - removing intermediate files..."
    rm -f "${VCF_CALL}.vcf.gz" "${VCF_CALL}.vcf.gz.csi"
else
    echo "H.calls file is empty or missing! Skipping cleanup."
fi

echo "Done chr ${CHRNAME}."
