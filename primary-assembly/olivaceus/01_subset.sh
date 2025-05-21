#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name Polv_sort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH -o Polv_sort.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

# Sort the HiFi reads in ascending order
seqkit sort --by-length --reverse P_olivaceus_HiFi_reads.fasta > sorted_hifi.fasta

# Take only 30% of the longest reads
seqkit head -n 1800000 sorted_hifi.fasta -o P_olivaceus_HiFi_top30.fasta

# went from 78G to 36G
