#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --partition=batch_long
#SBATCH --job-name braker
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=120:00:00
#SBATCH -o braker_PcSW.%j.out
#SBATCH -A lp_svbelleghem

#BRAKER and AUGUSTUS are both installed on the VSC_DATA folder
#BRAKER was installed with the singularity build
#AUGUSTUS was installed using Docker
#AUGUSTUS folder is moved INTO the BRAKER folder (that way BRAKER can access the config file)

export BRAKER_SIF=/data/leuven/357/vsc35707/BRAKER/braker3.sif

singularity exec \
  -B /scratch/leuven/357/vsc35707/annotation \
  -B /data/leuven/357/vsc35707/BRAKER/Augustus/config \
  "$BRAKER_SIF" \
  braker.pl \
  --threads=24 \
  --AUGUSTUS_CONFIG_PATH=/data/leuven/357/vsc35707/BRAKER/Augustus/config \
  --species=P_chalceus \
  --genome=Pchalceus_SW.sorted.fasta.masked \
  --prot_seq=coleoptera_proteins.clean.fasta \
  --bam=merged_bams/all_rnaseq_merged.bam \
  --softmasking
