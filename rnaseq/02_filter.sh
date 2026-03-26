#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --partition=batch_icelake
#SBATCH --job-name filter_rna
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o filter_rna.%A_%a.out
#SBATCH --array=1-12

THREADS=12

# Input and output directories
BAM_DIR="hisat2_bams"
OUT_DIR="filtered_bams"

mkdir -p "$OUT_DIR"

# Activate environment
conda activate variant_tools

# -------------------------
# Sample list
# -------------------------
samples=(LW_01 LW_02 LW_03 LW_04 LW_05 LW_06 \
SW_01 SW_02 SW_03 SW_04 SW_05 SW_06)

ID=$((SLURM_ARRAY_TASK_ID - 1))
sample="${samples[$ID]}"

INPUT_BAM="${BAM_DIR}/${sample}.sorted.bam"
FILTERED_BAM="${OUT_DIR}/${sample}.filtered.bam"
FINAL_BAM="${OUT_DIR}/${sample}.filtered.sorted.bam"

echo "=============================="
echo "Processing sample: $sample"
echo "Input BAM: $INPUT_BAM"
echo "=============================="

# Check file exists
if [[ ! -s "$INPUT_BAM" ]]; then
    echo "ERROR: missing BAM $INPUT_BAM"
    exit 1
fi

# -------------------------
# Filter reads
# -------------------------
echo "Filtering reads..."

samtools view \
  -@ "$THREADS" \
  -f 0x02 \
  -q 30 \
  -b "$INPUT_BAM" \
  > "$FILTERED_BAM"

# -------------------------
# Sort BAM
# -------------------------
echo "Sorting BAM..."

samtools sort \
  -@ "$THREADS" \
  "$FILTERED_BAM" \
  -o "$FINAL_BAM"

# -------------------------
# Index BAM
# -------------------------
echo "Indexing BAM..."

samtools index "$FINAL_BAM"

# Remove intermediate file
rm "$FILTERED_BAM"

echo "Finished processing $sample"
