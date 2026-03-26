#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --partition=batch_icelake
#SBATCH --job-name hisat2_rna
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o hisat2_rna.%A_%a.out
#SBATCH --array=1-12

# -------------------------
# settings
# -------------------------
THREADS=12
REF="Pchalceus_SW.sorted.fasta"
INDEX_PREFIX="Pchalceus_SW_hisat2"
READS_DIR="/scratch/leuven/357/vsc35707/rna-seq/reads"
OUT_DIR="hisat2_bams"

mkdir -p "$OUT_DIR"

# Activate environment 
conda activate variant_tools

# -------------------------
# Sample list
# -------------------------
samples=(LW_01 LW_02 LW_03 LW_04 LW_05 LW_06 SW_01 SW_02 SW_03 SW_04 SW_05 SW_06)

ID=$((SLURM_ARRAY_TASK_ID - 1))
sample="${samples[$ID]}"

R1="${READS_DIR}/${sample}_R1.fq.gz"
R2="${READS_DIR}/${sample}_R2.fq.gz"

echo "=============================="
echo "Sample: $sample"
echo "R1: $R1"
echo "R2: $R2"
echo "=============================="

# -------------------------
# Check input files
# -------------------------
if [[ ! -s "$R1" ]]; then
    echo "ERROR: missing file $R1" >&2
    exit 1
fi

if [[ ! -s "$R2" ]]; then
    echo "ERROR: missing file $R2" >&2
    exit 1
fi

# -------------------------
# Make sure HISAT2 index is built
# -------------------------

# -------------------------
# Map reads
# -------------------------
echo "Mapping $sample with HISAT2..."

hisat2 \
    -p "$THREADS" \
    -x "$INDEX_PREFIX" \
    -1 "$R1" \
    -2 "$R2" \
    --summary-file "${OUT_DIR}/${sample}.hisat2.summary.txt" \
| samtools view -bS - \
| samtools sort -@ "$THREADS" -o "${OUT_DIR}/${sample}.sorted.bam" -

# -------------------------
# Index BAM
# -------------------------
samtools index "${OUT_DIR}/${sample}.sorted.bam"

echo "Finished $sample"
