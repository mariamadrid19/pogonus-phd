#!/bin/bash -l 
#SBATCH --cluster=genius
#SBATCH --job-name call_snps
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o call_snps.%j.out
#SBATCH --array=1-11

# Index starts at 0
ID=$((SLURM_ARRAY_TASK_ID - 1))

# Load required modules
module load BCFtools/1.12-GCC-10.3.0 
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate variant_tools

# Reference files
REFNAME="Pchal_Bar_SW"
REF="/scratch/leuven/357/vsc35707/popgen/new-genome-july-2025/Pchalceus_SW.sorted.fasta"
FAI="${REF}.fai"

# Output directories
VCF_DIR="vcfs"
OUT_DIR="final-calls"
mkdir -p "$VCF_DIR"
mkdir -p "$OUT_DIR"

# Sample IDs (all samples, 103)
samples=(\
GC129388 GC129389 GC129390 GC129391 GC129392 \
GC129393 GC129394 GC129395 GC129396 GC129397 \
GC129398 GC129399 GC129400 GC129401 GC129402 \
GC129403 GC129404 GC129405 GC129406 GC129407 \
GC129408 GC129409 GC129410 GC129411 GC129412 \
GC129413 GC129414 GC129415 GC129416 GC129417 \
GC129418 GC129419 GC129420 GC129421 GC129422 \
GC129423 GC129424 GC129425 GC129426 GC129427 \
GC129428 GC129429 GC129430 GC129431 GC129432 \
GC129433 GC129434 GC129435 GC129437 GC129438 \
GC129439 GC129440 GC136078 GC136079 GC136080 \
GC136081 GC136082 GC136083 GC136084 GC136085 \
GC136086 GC136087 GC136088 GC136089 GC136090 \
GC136091 GC136092 GC136093 GC136094 GC136095 \
GC136096 GC136097 GC136098 GC136099 GC136100 \
GC136101 GC136102 GC136103 GC136104 GC136105 \
GC136106 GC136107 GC136108 GC136109 GC136110 \
GC136111 GC136112 GC136113 GC136114 GC136115 \
GC136116 GC136117 GC136118 GC136119 GC136120 \
GC136121 GC136122 GC136123 GC136124 GC136125 \
GC136126 GC136127 GC136128)

CHR=(CHR1 CHR2 CHR3 CHR4 CHR5 CHR6 CHR7 CHR8 CHR9 CHR10 CHR11)
names=(1 2 3 4 5 6 7 8 9 10 11)

# Build BAM file list from 'bams/' directory
ALL_LIST=""
for SAMPLE in "${samples[@]}"; do
    ALL_LIST+=" bams/${SAMPLE}.${REFNAME}.filtered.sorted.dedup.bam"
done

# Step 1: Call variants with FORMAT annotations
bcftools mpileup -Oz --threads 36 --fasta-ref "$REF" --regions $(echo "${CHR[ID]}") $ALL_LIST --annotate FORMAT/DP | bcftools call --multiallelic-caller -Oz -f GQ -o $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").vcf.gz

# Step 2: Filter with vcftools
vcftools --gzvcf $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").vcf.gz --recode --remove-indels --stdout | bgzip > $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").filtered.vcf.gz

# Step 4: Parse VCF with custom script
/data/leuven/357/vsc35707/miniconda3/bin/python3.11 parseVCF.py \
  --gtf flag=GQ min=30 gtTypes=Het \
  --gtf flag=GQ min=30 gtTypes=HomAlt \
  --gtf flag=DP min=10 \
  --skipIndels \
  -i $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").filtered.vcf.gz \
  | gzip > $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").calls.gz

echo "Calls file created: $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").calls.gz" 

zcat $VCF_DIR/$REFNAME.chr_$(echo "${names[ID]}").calls.gz | awk -v OFS="\t" 'NR==1 {
    for (i = 1; i <= NF; i++) {
        gsub(/^bams\//, "", $i);
        gsub(/\.Pchal_Bar_SW\.filtered\.sorted\.dedup\.bam$/, "", $i);
    }
}
{ print }
' | gzip -c > $OUT_DIR/$REFNAME.chr_$(echo "${names[ID]}").H.calls.gz

echo "Header fixed and saved: $OUT_DIR/$REFNAME.chr_$(echo "${names[ID]}").H.calls.gz"
