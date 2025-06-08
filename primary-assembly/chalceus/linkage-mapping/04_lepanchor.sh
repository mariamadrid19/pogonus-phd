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

cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap
