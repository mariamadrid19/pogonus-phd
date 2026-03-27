#!/bin/bash -l
#SBATCH --job-name=eggnog
#SBATCH --cluster=genius
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=48:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/03_eggnog_%j.out
#SBATCH -e /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/03_eggnog_%j.err

set -euo pipefail

PROJECT=/scratch/leuven/357/vsc35707/annotation/func-annotation
INPUT=$PROJECT/results/01_longest_isoforms/braker.longest_isoforms.aa
OUTDIR=$PROJECT/results/03_eggnog
DBDIR=$PROJECT/databases/eggnog
TMPDIR_LOCAL=$PROJECT/tmp/eggnog_${SLURM_JOB_ID}

mkdir -p "$OUTDIR"
mkdir -p "$TMPDIR_LOCAL"

echo "======================================"
echo "eggNOG-mapper run"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Input: $INPUT"
echo "Output dir: $OUTDIR"
echo "DB dir: $DBDIR"
echo "Temp dir: $TMPDIR_LOCAL"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "======================================"

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate eggnog

echo ""
echo "Versions:"
which emapper.py
which diamond
emapper.py -v || true
diamond --version || true

echo ""
echo "Checking required files..."
[[ -s "$INPUT" ]] || { echo "ERROR: input file not found: $INPUT"; exit 1; }
[[ -s "$DBDIR/eggnog.db" ]] || { echo "ERROR: missing $DBDIR/eggnog.db"; exit 1; }
[[ -s "$DBDIR/eggnog.taxa.db" ]] || { echo "ERROR: missing $DBDIR/eggnog.taxa.db"; exit 1; }
[[ -s "$DBDIR/eggnog_proteins.dmnd" ]] || { echo "ERROR: missing $DBDIR/eggnog_proteins.dmnd"; exit 1; }

echo ""
echo "Starting eggNOG-mapper..."
emapper.py \
  -i "$INPUT" \
  --itype proteins \
  -m diamond \
  --tax_scope eukaryota \
  --cpu "${SLURM_CPUS_PER_TASK}" \
  --data_dir "$DBDIR" \
  --temp_dir "$TMPDIR_LOCAL" \
  --output braker_longest_isoforms \
  --output_dir "$OUTDIR"

echo ""
echo "eggNOG-mapper finished: $(date)"
echo ""
echo "Output files:"
ls -lh "$OUTDIR"

echo ""
echo "Cleaning temp directory..."
rm -rf "$TMPDIR_LOCAL"

echo ""
echo "Done."
