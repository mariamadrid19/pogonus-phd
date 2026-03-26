#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=star_rna
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o star_rna.%A_%a.out
#SBATCH --array=0-14

conda activate thesis

READS_DIR="/scratch/leuven/357/vsc35707/annotation/reads"
STAR_INDEX="/scratch/leuven/357/vsc35707/annotation/STAR_index_softmasked"
OUTDIR="/scratch/leuven/357/vsc35707/annotation/star_bams"
THREADS=24

mkdir -p "$OUTDIR"

samples=(
  "LW_01" "LW_02" "LW_03" "LW_04" "LW_05" "LW_06"
  "SW_01" "SW_02" "SW_03" "SW_04" "SW_05" "SW_06"
  "SRR424340" "SRR424342" "SRR424344"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

if [[ "$sample" == SRR* ]]; then
    R1="${READS_DIR}/${sample}_1.fastq.gz"
    R2="${READS_DIR}/${sample}_2.fastq.gz"
else
    R1="${READS_DIR}/${sample}_R1.fq.gz"
    R2="${READS_DIR}/${sample}_R2.fq.gz"
fi

if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: missing reads for $sample"
    echo "R1=$R1"
    echo "R2=$R2"
    exit 1
fi

PREFIX="${OUTDIR}/${sample}."

STAR \
  --runThreadN "$THREADS" \
  --genomeDir "$STAR_INDEX" \
  --readFilesIn "$R1" "$R2" \
  --readFilesCommand zcat \
  --outFileNamePrefix "$PREFIX" \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic

samtools index "${PREFIX}Aligned.sortedByCoord.out.bam"

echo "Finished STAR mapping for $sample"
