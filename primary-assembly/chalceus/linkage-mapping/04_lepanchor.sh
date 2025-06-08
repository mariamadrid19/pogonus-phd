#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name lepanchor
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o lepanchor.%j.out

# Change to the working directory
cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap

# Define directories
LEPANCHOR="/data/leuven/357/vsc35707/LepAnchor"

# Clean map file
java -cp bin/ CleanMap map=map.txt >map.clean
