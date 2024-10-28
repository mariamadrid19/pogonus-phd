#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name braker
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o braker.%j.out
#SBATCH -A lp_svbelleghem

#BRAKER and AUGUSTUS are both installed in the PWD (locally, on the SCRATCH)
#BRAKER was installed with the make file
#AUGUSTUS was installed using Docker (as job, requires a lot of memory)

singularity exec -B ${PWD}:${PWD} $PWD/BRAKER/braker3.sif braker.pl --threads=24 --AUGUSTUS_CONFIG_PATH=$PWD/Augustus/config --species=P_chalceus --genome=sorted_prim_dud.fasta.masked --bam=POG_mapped_RNA_dud.sorted.filtered.bam --softmasking
