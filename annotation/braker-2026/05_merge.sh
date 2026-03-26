#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=merge_rna
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=08:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o merge_rna.%j.out

conda activate thesis

STAR_DIR="/scratch/leuven/357/vsc35707/annotation/star_bams"
ISOSEQ_BAM="/scratch/leuven/357/vsc35707/annotation/isoseq_bams/POG_IsoSeq_HiFi.sorted.bam"
OUTDIR="/scratch/leuven/357/vsc35707/annotation/merged_bams"
THREADS=20

mkdir -p "$OUTDIR"

# Check inputs exist
ls "${STAR_DIR}"/*.Aligned.sortedByCoord.out.bam > /dev/null
[[ -s "$ISOSEQ_BAM" ]]

samtools merge -@ "$THREADS" \
  "${OUTDIR}/all_rnaseq_merged.bam" \
  "${STAR_DIR}"/*.Aligned.sortedByCoord.out.bam \
  "$ISOSEQ_BAM"

samtools index "${OUTDIR}/all_rnaseq_merged.bam"
