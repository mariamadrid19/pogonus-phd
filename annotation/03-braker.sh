#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name braker
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=48:00:00
#SBATCH -o braker.%j.out
#SBATCH -A lp_svbelleghem

conda activate braker

braker.pl --cores=12 --species=P_chalceus --genome=sorted_prim_dud.soft_masked.fasta --bam=POG_mapped_RNA_dud.sorted.filtered.bam --softmasking --AUGUSTUS_CONFIG_PATH=/data/leuven/357/vsc35707/Augustus/config/
