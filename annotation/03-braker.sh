#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name braker
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o braker.%j.out
#SBATCH -A lp_svbelleghem

#BRAKER and AUGUSTUS are both installed on the VSC_DATA folder
#BRAKER was installed with the singularity build
#AUGUSTUS was installed using Docker
#Both were sent as jobs, installation requires a lot of memory 
#AUGUSTUS folder is moved INTO the BRAKER folder (that way BRAKER can access the config file)

#to execute braker3.sif
export BRAKER_SIF=/data/leuven/357/vsc35707/BRAKER/braker3.sif

singularity exec -B /data/leuven/357/vsc35707/BRAKER/braker3.sif braker.pl --threads=24 --AUGUSTUS_CONFIG_PATH=/data/leuven/357/vsc35707/BRAKER/Augustus/config --species=P_chalceus --genome=sorted_prim_dud.fasta.masked --bam=POG_mapped_RNA_dud.sorted.filtered.bam --softmasking
