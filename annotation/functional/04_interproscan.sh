#!/bin/bash -l
#SBATCH --job-name=interpro
#SBATCH --cluster=wice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/04_interpro_%j.out
#SBATCH -e /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/04_interpro_%j.err

set -euo pipefail

module load Java/21.0.8

PROJECT=/scratch/leuven/357/vsc35707/annotation/func-annotation
INPUT_RAW=$PROJECT/results/01_longest_isoforms/braker.longest_isoforms.aa
INPUT_CLEAN=$PROJECT/results/01_longest_isoforms/braker.longest_isoforms.no_stop.aa
OUTDIR=$PROJECT/results/04_interproscan
TMPDIR_LOCAL=$PROJECT/tmp/interproscan_${SLURM_JOB_ID}
INTERPRO=/scratch/leuven/357/vsc35707/my_interproscan/interproscan-5.77-108.0/interproscan.sh

mkdir -p "$OUTDIR"
mkdir -p "$TMPDIR_LOCAL"

echo "======================================"
echo "InterProScan run"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Raw input: $INPUT_RAW"
echo "Clean input: $INPUT_CLEAN"
echo "Output dir: $OUTDIR"
echo "Temp dir: $TMPDIR_LOCAL"
echo "InterProScan: $INTERPRO"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "======================================"

echo ""
echo "Checking required files..."
[[ -s "$INPUT_RAW" ]] || { echo "ERROR: input file not found: $INPUT_RAW"; exit 1; }
[[ -x "$INTERPRO" ]] || { echo "ERROR: InterProScan not executable: $INTERPRO"; exit 1; }

echo ""
echo "Creating cleaned FASTA without '*' characters..."
awk '
BEGIN{FS=OFS=""}
(/^>/){print; next}
{gsub(/\*/,""); print}
' "$INPUT_RAW" > "$INPUT_CLEAN"

echo ""
echo "Checking for remaining '*' characters..."
if grep -q '\*' "$INPUT_CLEAN"; then
    echo "ERROR: cleaned FASTA still contains '*' characters"
    exit 1
fi

echo ""
echo "Starting InterProScan..."
"$INTERPRO" \
  -i "$INPUT_CLEAN" \
  -f tsv,gff3 \
  -goterms \
  -iprlookup \
  -pa \
  -cpu "${SLURM_CPUS_PER_TASK}" \
  --tempdir "$TMPDIR_LOCAL" \
  -b "$OUTDIR/braker_longest_isoforms"

echo ""
echo "InterProScan finished: $(date)"
echo ""
echo "Output files:"
ls -lh "$OUTDIR"

echo ""
echo "Cleaning temp directory..."
rm -rf "$TMPDIR_LOCAL"

echo ""
echo "Done."
