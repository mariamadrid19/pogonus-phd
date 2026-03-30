#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name merge_annotations
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=2:00:00
#SBATCH -o /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/merge.%j.out
#SBATCH --error=/scratch/leuven/357/vsc35707/annotation/func-annotation/logs/merge.%j.err
#SBATCH -A lp_svbelleghem

module load Python/3.13.1-GCCcore-14.2.0

python3 /scratch/leuven/357/vsc35707/annotation/func-annotation/scripts/merge_annotations.py
