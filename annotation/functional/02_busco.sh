#!/bin/bash -l
#SBATCH --job-name=busco_proteins
#SBATCH --cluster=genius
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=2:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/02_busco_%j.out
#SBATCH -e /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/02_busco_%j.err

set -euo pipefail

PROJECT=/scratch/leuven/357/vsc35707/annotation/func-annotation
INPUT=$PROJECT/results/01_longest_isoforms/braker.longest_isoforms.aa
OUTDIR=$PROJECT/results/02_busco
LINEAGE=insecta_odb12

mkdir -p "$OUTDIR"

echo "=============================="
echo "BUSCO run"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Input: $INPUT"
echo "Output dir: $OUTDIR"
echo "Lineage: $LINEAGE"
echo "=============================="

conda activate busco

busco \
  -i "$INPUT" \
  -m proteins \
  -l "$LINEAGE" \
  -o braker_longest_isoforms_busco \
  --out_path "$OUTDIR" \
  -c ${SLURM_CPUS_PER_TASK}

echo "=============================="
echo "BUSCO finished"
echo "Date: $(date)"
echo "=============================="
