#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name map
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --time=60:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o map_dud.%A_%a.out
#SBATCH --array=1-185

# Map + filter + sort + index for all samples listed in samples.present.txt

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID - 1))

conda activate variant_tools
# samtools, bwa are here
# make sure to index the reference fasta beforehand with bwa index

echo "================="
echo "SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"

# ---- Sample list (only samples that actually exist) ----
SAMPLES_LIST="samples.present.txt"
if [[ ! -f "$SAMPLES_LIST" ]]; then
  echo "ERROR: Cannot find $SAMPLES_LIST in $(pwd)" >&2
  echo "Submit from the reads/ directory or set SAMPLES_LIST to an absolute path." >&2
  exit 1
fi

mapfile -t samples < "$SAMPLES_LIST"

# Safety check: array range vs list length
if [[ "$SLURM_ARRAY_TASK_ID" -gt "${#samples[@]}" ]]; then
  echo "ERROR: Task ID ${SLURM_ARRAY_TASK_ID} exceeds number of samples (${#samples[@]})" >&2
  exit 1
fi

echo "${samples[ID]}"

# ---- Some folder and file paths to use later ----
REF=/scratch/leuven/357/vsc35707/GWAS/Dudzele/Pchalceus_SW.sorted.fasta
REFNAME=Pchal_Bar_SW
BWAout=/scratch/leuven/357/vsc35707/GWAS/Dudzele/bams
READS_DIR=/scratch/leuven/357/vsc35707/GWAS/Dudzele/reads

mkdir -p "$BWAout"

FILE1="$READS_DIR/${samples[ID]}_R1.fq.gz"
FILE2="$READS_DIR/${samples[ID]}_R2.fq.gz"

# Fail fast if FASTQs missing
if [[ ! -s "$FILE1" || ! -s "$FILE2" ]]; then
  echo "ERROR: Missing FASTQ(s) for ${samples[ID]}" >&2
  echo "  FILE1=$FILE1" >&2
  echo "  FILE2=$FILE2" >&2
  exit 1
fi

echo "Mapping reads: ${samples[ID]}"
# Map reads using bwa mem
bwa mem -t 36 -M "$REF" "$FILE1" "$FILE2" \
  | samtools view -bS - \
  > "$BWAout/${samples[ID]}.$REFNAME.bam"

echo "Filtering reads: ${samples[ID]}"
# Filter using samtools
samtools view -f 0x02 -q 20 -b \
  "$BWAout/${samples[ID]}.$REFNAME.bam" \
  > "$BWAout/${samples[ID]}.$REFNAME.filtered.bam"

# Remove .bam if .filtered.bam is larger than 500 KB
if [[ -f "$BWAout/${samples[ID]}.$REFNAME.filtered.bam" ]]; then
  if [[ "$(stat -c %s "$BWAout/${samples[ID]}.$REFNAME.filtered.bam")" -gt 512000 ]]; then
    rm -f "$BWAout/${samples[ID]}.$REFNAME.bam"
  fi
fi

echo "Sorting reads: ${samples[ID]}"
# Sort using samtools
samtools sort \
  "$BWAout/${samples[ID]}.$REFNAME.filtered.bam" \
  -o "$BWAout/${samples[ID]}.$REFNAME.filtered.sorted.bam"

# Remove .filtered.bam if .filtered.sorted.bam is larger than 500 KB
if [[ -f "$BWAout/${samples[ID]}.$REFNAME.filtered.sorted.bam" ]]; then
  if [[ "$(stat -c %s "$BWAout/${samples[ID]}.$REFNAME.filtered.sorted.bam")" -gt 512000 ]]; then
    rm -f "$BWAout/${samples[ID]}.$REFNAME.filtered.bam"
  fi
fi

samtools index "$BWAout/${samples[ID]}.$REFNAME.filtered.sorted.bam"

echo "Finished with ${samples[ID]}"
