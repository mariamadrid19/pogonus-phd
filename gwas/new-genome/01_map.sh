#!/bin/bash -l 
#SBATCH --cluster=genius
#SBATCH --job-name map
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=36
#SBATCH --time=60:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o map_all.%j.out
#SBATCH --array=1-192

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

conda activate variant_tools
# samtools, bwa, vcftools are here 

echo "================="

# Sample IDs (all samples, 192)
samples=()
for i in $(seq -w 1 192); do
  samples+=("Pc25Np${i}")
done

echo "${samples[ID]}"

# Some folder and file paths to use later
REF=/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/Pchalceus_SW.sorted.fasta
REFNAME=Pchal_Bar_SW
BWAout=/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/bams
FILE1=/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/reads/$(echo "${samples[ID]}")_R1.fq.gz
FILE2=/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/reads/$(echo "${samples[ID]}")_R2.fq.gz

echo "Mapping reads: ${samples[ID]}"
# Map reads using bwa mem
bwa mem -t 36 -M $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[ID]}").$REFNAME.bam

echo "Filtering reads: ${samples[ID]}"
# Filter using samtools
samtools view -f 0x02 -q 20 -b $BWAout/$(echo "${samples[ID]}").$REFNAME.bam > $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam

# Remove .bam if .filtered.bam is larger than 500 KB
if [ $(stat -c %s "$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam") -gt 512000 ]; then
    rm "$BWAout/$(echo "${samples[ID]}").$REFNAME.bam"
fi

echo "Sorting reads: ${samples[ID]}"
# Sort using samtools
samtools sort $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam -o $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

# Remove .filtered.bam if .filtered.sorted.bam is larger than 500 KB
if [ $(stat -c %s "$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam") -gt 512000 ]; then
    rm "$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam"
fi

samtools index $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

echo "Finished with ${samples[ID]}"
