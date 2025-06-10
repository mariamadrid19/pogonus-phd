#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name call
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o call.%j.out

# Load modules
module load SAMtools/1.16.1-GCC-11.3.0
module load Java/21.0.2 

# Change to the working directory
cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap

# Define directories
LEPMAP="/data/leuven/357/vsc35707/LepMap3"

# Run samtools mpileup and calculate posterior genotypes
samtools mpileup -q 20 -Q 10 -s $(cat sorted_bams.txt) \
| java -cp $LEPMAP/bin/ Pileup2Likelihoods mappingFile=mapping.txt minCoverage=10 numLowerCoverage=5 | gzip > post.gz
