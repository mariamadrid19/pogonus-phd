#!/bin/bash -l
#SBATCH --job-name=diamond
#SBATCH --cluster=wice
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=36:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/05_diamond_%j.out
#SBATCH -e /scratch/leuven/357/vsc35707/annotation/func-annotation/logs/05_diamond_%j.err

set -euo pipefail

PROJECT=/scratch/leuven/357/vsc35707/annotation/func-annotation
INPUT=$PROJECT/results/01_longest_isoforms/braker.longest_isoforms.no_stop.aa
OUTDIR=$PROJECT/results/05_diamond
DBDIR=$PROJECT/databases/uniprot_swissprot
DBFASTA=$DBDIR/uniprot_sprot.fasta
DBNAME=$DBDIR/uniprot_sprot

mkdir -p "$OUTDIR"

echo "======================================"
echo "DIAMOND vs Swiss-Prot"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "Input: $INPUT"
echo "Output dir: $OUTDIR"
echo "DB fasta: $DBFASTA"
echo "DB base: $DBNAME"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "======================================"

echo ""
echo "Checking required files..."
[[ -s "$INPUT" ]] || { echo "ERROR: input file not found: $INPUT"; exit 1; }
[[ -s "$DBFASTA" ]] || { echo "ERROR: Swiss-Prot FASTA not found: $DBFASTA"; exit 1; }

echo ""
echo "Loading DIAMOND..."
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate eggnog

which diamond
diamond --version

echo ""
echo "Building DIAMOND database if needed..."
if [[ ! -s "${DBNAME}.dmnd" ]]; then
    diamond makedb \
        --in "$DBFASTA" \
        -d "$DBNAME"
fi

echo ""
echo "Running DIAMOND blastp..."
diamond blastp \
    --db "${DBNAME}.dmnd" \
    --query "$INPUT" \
    --out "$OUTDIR/braker_longest_isoforms_vs_swissprot.tsv" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
    --max-target-seqs 5 \
    --evalue 1e-5 \
    --threads "${SLURM_CPUS_PER_TASK}"

echo ""
echo "Finished: $(date)"
echo ""
echo "Output files:"
ls -lh "$OUTDIR"

echo ""
echo "Done."
