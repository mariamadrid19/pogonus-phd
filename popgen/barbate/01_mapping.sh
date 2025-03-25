#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name barbate_map
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o barbate_map.%j.out
#SBATCH --array=1-21

cd /scratch/leuven/357/vsc35707/BAR_mapping

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BWA/0.7.17-foss-2018a
module load SAMtools/1.18-GCC-12.3.0
module load picard/2.18.23-Java-1.8.0_171
module load minimap2/2.26-GCCcore-12.3.0
module load Python/3.7.0-foss-2018a
module load tabix/0.2.6-GCCcore-6.4.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

echo "================="

# Sample IDs (20 samples)
samples=(Bar2_01 Bar2_02 Bar2_03 Bar2_04 Bar2_05 Bar2_06 Bar2_07 Bar2_08 Bar2_09 Bar2_10\
	Bar4_01 Bar4_02 Bar4_03 Bar4_04 Bar4_05 Bar4_06 Bar4_07 Bar4_08 Bar4_09 Bar4_10)

echo "${samples[ID]}"

# Some folder and file paths to use later
REF=/scratch/leuven/357/vsc35707/BAR_mapping/sorted_prim_dud.fasta
REFNAME=dudPrim
BWAout=/scratch/leuven/357/vsc35707/BAR_mapping/bams
FILE1=/scratch/leuven/357/vsc35707/BAR_mapping/$(echo "${samples[ID]}")_R1.fastq.gz
FILE2=/scratch/leuven/357/vsc35707/BAR_mapping/$(echo "${samples[ID]}")_R2.fastq.gz

# Map reads using bwa mem
bwa mem -t 20 -M $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[ID]}").$REFNAME.bam

# Filter using samtools
samtools view -f 0x02 -q 20 -b $BWAout/$(echo "${samples[ID]}").$REFNAME.bam > $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam

# Sort using samtools
samtools sort $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam -o $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

# Remove PCR duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
-INPUT $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam \
-OUTPUT $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam \
-REMOVE_DUPLICATES true \
-METRICS_FILE $BWAout/$(echo "${samples[ID]}").$REFNAME.dup_metrics.txt \
-ASSUME_SORTED true

# Remove intermediate files
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam 

samtools index $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam
