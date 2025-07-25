#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name Spain_map
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o Spain_map.%j.out
#SBATCH --array=1-10

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BWA/0.7.17-GCC-10.3.0
module load SAMtools/1.16.1-GCC-11.3.0

echo "================="

# Sample IDs (10 samples)
samples=(GC136107 GC136108 GC136109 GC136110 GC136111 GC136112 GC136113 GC136114 GC136115 GC136116)

echo "${samples[ID]}"

# Some folder and file paths to use later
REF=/scratch/leuven/357/vsc35707/popgen/P_chalceus_REF1.fa
REFNAME=REF1
BWAout=/scratch/leuven/357/vsc35707/popgen/bams
FILE1=/scratch/leuven/357/vsc35707/popgen/$(echo "${samples[ID]}")_R1.fastq.gz
FILE2=/scratch/leuven/357/vsc35707/popgen/$(echo "${samples[ID]}")_R2.fastq.gz

# Check the directory exists
mkdir -p $BWAout

# Map reads using bwa mem
bwa mem -t 20 -M $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[ID]}").$REFNAME.bam

# Filter using samtools
samtools view -f 0x02 -q 20 -b $BWAout/$(echo "${samples[ID]}").$REFNAME.bam > $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam

# Sort using samtools
samtools sort $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam -o $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

module load picard/2.18.23-Java-1.8.0_171

# Remove PCR duplicates
java -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam \
-O $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.dedup.bam \
-REMOVE_DUPLICATES true -M $BWAout/$(echo "${samples[ID]}").$REFNAME.dup_metrics.txt -ASSUME_SORTED true

# Remove intermediate files
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam 

samtools index $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.dedup.bam
