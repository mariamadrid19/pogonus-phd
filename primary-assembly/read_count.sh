#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name reads 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:10:00 
#SBATCH -A lp_svbelleghem

# Define the directory containing FASTQ files
FASTQ_DIR=/scratch/leuven/357/vsc35707

# Loop through all .fastq.gz files in the directory
for file in "$FASTQ_DIR"/*.fastq.gz; do
    if [[ -f "$file" ]]; then
        # Count the total number of lines and divide by 4 (each read spans 4 lines)
        READ_COUNT=$(zcat "$file" | wc -l | awk '{print $1 / 4}')
        echo "$file: $READ_COUNT reads"
    fi
done 
