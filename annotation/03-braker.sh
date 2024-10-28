#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name braker
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o braker.%j.out
#SBATCH -A lp_svbelleghem

conda activate braker

export AUGUSTUS_CONFIG_PATH=/data/leuven/357/vsc35707/Augustus/config

braker3.sif braker.pl --threads=24 --species=P_chalceus --genome=sorted_prim_dud.fasta.masked --bam=POG_mapped_RNA_dud.sorted.filtered.bam --softmasking
