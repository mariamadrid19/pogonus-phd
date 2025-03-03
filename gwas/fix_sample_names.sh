#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name fix 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=20:00:00 
#SBATCH -A lp_svbelleghem

# Create output directory if it doesn't exist
mkdir -p fixed

# Process R1 files
for file in *_R1.fastq.gz; do
    zcat "$file" | sed 's/^\(@.*\)_1/\1/' | gzip > "fixed/$file"
    echo "Processed R1: $file -> fixed/$file"
done

# Process R2 files
for file in *_R2.fastq.gz; do
    zcat "$file" | sed 's/^\(@.*\)_2/\1/' | gzip > "fixed/$file"
    echo "Processed R2: $file -> fixed/$file"
done

echo "All files processed successfully! Fixed files are in the 'fixed/' directory."
