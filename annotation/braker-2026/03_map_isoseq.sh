#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=isoseq_map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o isoseq_map.%j.out

conda activate thesis

GENOME="/scratch/leuven/357/vsc35707/annotation/Pchalceus_SW.sorted.fasta.masked"
ISOSEQ="/scratch/leuven/357/vsc35707/annotation/reads/POG_IsoSeq_HiFi_demux.fastq"
OUTDIR="/scratch/leuven/357/vsc35707/annotation/isoseq_bams"
THREADS=24

mkdir -p "$OUTDIR"

OUTBAM="${OUTDIR}/POG_IsoSeq_HiFi.sorted.bam"

minimap2 \
  -t "$THREADS" \
  -ax splice:hq \
  -uf \
  --secondary=no \
  "$GENOME" "$ISOSEQ" \
  | samtools sort -@ "$THREADS" -o "$OUTBAM"

samtools index "$OUTBAM"
