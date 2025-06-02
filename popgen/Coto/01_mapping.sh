#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=20:00:00 
#SBATCH -A lp_svbelleghem

#Script written by Steven Van Belleghem (2024)

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BWA
module load SAMtools
