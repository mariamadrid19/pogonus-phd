#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=featureCounts
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o featureCounts.%j.out
#SBATCH -e featureCounts.%j.err

set -euo pipefail

THREADS=12

BAM_DIR="/scratch/leuven/357/vsc35707/rna-seq/filtered_bams"
OUT_DIR="/scratch/leuven/357/vsc35707/rna-seq/counts"
GTF="/scratch/leuven/357/vsc35707/annotation/func-annotation/results/08_minimal_annotation/braker.minimal.annotated.gtf"

mkdir -p "$OUT_DIR"

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate variant_tools

featureCounts \
  -T "$THREADS" \
  -p \
  --countReadPairs \
  -B \
  -C \
  -t exon \
  -g gene_id \
  -a "$GTF" \
  -o "$OUT_DIR/gene_counts.txt" \
  "$BAM_DIR"/*.filtered.sorted.bam
