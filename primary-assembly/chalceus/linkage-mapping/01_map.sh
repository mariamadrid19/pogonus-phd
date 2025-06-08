#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa_map
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o bwa_map.%j.out
#SBATCH --array=1-80

# Load modules
module load BWA/0.7.17-GCC-10.3.0
module load SAMtools/1.16.1-GCC-11.3.0

# Change to the working directory
cd /scratch/leuven/357/vsc35707/linkage-mapping/

# Read sample name from file based on array task ID
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ./samples/samples.txt)

# Define genome file
GENOME="./genome/P_chalceus_broken.fa"

# Define input and output files
READ1="./reads/${SAMPLE}.1.fil.fq_1.gz"
READ2="./reads/${SAMPLE}.2.fil.fq_2.gz"
BAM_OUT="./bam/${SAMPLE}.bam"

mkdir -p ./bam

# Run BWA MEM, convert to BAM and sort
bwa mem -t 20 $GENOME $READ1 $READ2 | samtools view -bS | samtools sort -o $BAM_OUT
samtools index $BAM_OUT
