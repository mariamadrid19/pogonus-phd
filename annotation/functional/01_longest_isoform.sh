#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name isoform_extraction
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH -o /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/isoform_extraction.%j.out
#SBATCH --error=/scratch/leuven/357/vsc35707/annotation/func-annotation/logs/isoform_extraction.%j.err
#SBATCH -A lp_svbelleghem

module load Python/3.13.1-GCCcore-14.2.0

PROJECT=/scratch/leuven/357/vsc35707/annotation/func-annotation

python3 $PROJECT/scripts/extract_longest_isoforms.py \
  --gtf $PROJECT/input/braker/braker.gtf \
  --aa $PROJECT/input/braker/braker.aa \
  --cds $PROJECT/input/braker/braker.codingseq \
  --out-aa $PROJECT/results/01_longest_isoforms/braker.longest_isoforms.aa \
  --out-cds $PROJECT/results/01_longest_isoforms/braker.longest_isoforms.codingseq \
  --out-tsv $PROJECT/results/01_longest_isoforms/braker.longest_isoforms.tsv \
  2>&1 | tee $PROJECT/logs/01_extract_longest_isoforms.log
