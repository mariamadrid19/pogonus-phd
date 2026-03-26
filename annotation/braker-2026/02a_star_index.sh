#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=star_index
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o star_index.%j.out

conda activate thesis

GENOME="/scratch/leuven/357/vsc35707/annotation/Pchalceus_SW.sorted.fasta.masked"
STAR_INDEX="/scratch/leuven/357/vsc35707/annotation/STAR_index_softmasked"
THREADS=24

mkdir -p "$STAR_INDEX"

STAR \
  --runThreadN "$THREADS" \
  --runMode genomeGenerate \
  --genomeDir "$STAR_INDEX" \
  --genomeFastaFiles "$GENOME"
