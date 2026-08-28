#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=merge_annotations
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/leuven/357/vsc35707/annotation/func-annotation/logs/merge.%j.out
#SBATCH --error=/scratch/leuven/357/vsc35707/annotation/func-annotation/logs/merge.%j.err
#SBATCH --account=lp_svbelleghem

set -euo pipefail

module load Python/3.13.1-GCCcore-14.2.0

python3 /scratch/leuven/357/vsc35707/annotation/func-annotation/scripts/merge_annotations.py \
    /scratch/leuven/357/vsc35707/annotation/func-annotation
