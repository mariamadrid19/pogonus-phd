#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=map_contigs_SW
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=map_contigs_SW.%j.out
#SBATCH --error=map_contigs_SW.%j.err
#SBATCH --account=lp_svbelleghem

# ============================================================
# SW Pogonus Hi-C mapping pipeline
#
# Maps paired-end Arima Hi-C reads against the purge_dups
# purged SW primary assembly.
#
# Based on:
# https://github.com/ArimaGenomics/mapping_pipeline
# ============================================================

# ------------------------------------------------------------
# Main paths
# ------------------------------------------------------------

WORKDIR="/scratch/leuven/357/vsc35707/sw-assembly"
ASSEMBLY_DIR="${WORKDIR}/T2T_assembly"
HIC_DIR="${WORKDIR}/hiC_reads"
SCRIPT_DIR="${WORKDIR}/scripts"

REF="${ASSEMBLY_DIR}/purged.fa"
FAIDX="${REF}.fai"

HIC_R1="${HIC_DIR}/GC164167_TTCCAAGG-CCTTGTAG_S10_L001_R1_001.fastq.gz"
HIC_R2="${HIC_DIR}/GC164167_TTCCAAGG-CCTTGTAG_S10_L001_R2_001.fastq.gz"

# ------------------------------------------------------------
# Sample information
# ------------------------------------------------------------

SAMPLE="GC164167"
LABEL="Pogonus_chalceus_SW"

# Prefix used for the BWA index files
BWA_PREFIX="${ASSEMBLY_DIR}/purged"

REP_LABEL="${LABEL}_r"

MAPQ_FILTER=20
CPU="${SLURM_CPUS_PER_TASK:-36}"

# ------------------------------------------------------------
# Arima pipeline Perl scripts
# ------------------------------------------------------------

FILTER="${SCRIPT_DIR}/filter_five_end.pl"
COMBINER="${SCRIPT_DIR}/two_read_bam_combiner.pl"
STATS="${SCRIPT_DIR}/get_stats.pl"

# ------------------------------------------------------------
# Output directories
# ------------------------------------------------------------

MAP_DIR="${WORKDIR}/HiC_mapping_SW"

RAW_DIR="${MAP_DIR}/01_raw_bams"
FILT_DIR="${MAP_DIR}/02_filtered_bams"
TMP_DIR="${MAP_DIR}/03_temporary_files"
PAIR_DIR="${MAP_DIR}/04_paired_bams"
REP_DIR="${MAP_DIR}/05_deduplicated_bams"

mkdir -p \
    "$RAW_DIR" \
    "$FILT_DIR" \
    "$TMP_DIR" \
    "$PAIR_DIR" \
    "$REP_DIR"

# ------------------------------------------------------------
# Activate software environment
# ------------------------------------------------------------

source "/data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh"
conda activate thesis

# Avoid user-level Python packages interfering with Conda
unset PYTHONPATH
unset PYTHONHOME
export PYTHONNOUSERSITE=1

module purge
module load picard/3.4.0-Java-17
module load Java/21.0.8

# ------------------------------------------------------------
# Job information
# ------------------------------------------------------------

echo "============================================================"
echo "SW Hi-C mapping pipeline"
echo "============================================================"
echo "Node:             $(hostname)"
echo "Start time:       $(date)"
echo "Threads:          $CPU"
echo "Reference:        $REF"
echo "Hi-C R1:          $HIC_R1"
echo "Hi-C R2:          $HIC_R2"
echo "Output directory: $MAP_DIR"
echo

# ------------------------------------------------------------
# Check required files
# ------------------------------------------------------------

for file in \
    "$REF" \
    "$HIC_R1" \
    "$HIC_R2" \
    "$FILTER" \
    "$COMBINER" \
    "$STATS"
do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Required file is missing or empty:" >&2
        echo "$file" >&2
        exit 1
    fi
done

# Check required programs
for program in bwa samtools perl java
do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Program not found in PATH: $program" >&2
        exit 1
    fi
done

echo "Software:"
echo "BWA:      $(command -v bwa)"
echo "samtools: $(command -v samtools)"
echo "Perl:     $(command -v perl)"
echo "Picard:   ${EBROOTPICARD}/picard.jar"
echo

# ------------------------------------------------------------
# Step 0A: samtools FASTA index
# ------------------------------------------------------------

echo "============================================================"
echo "Step 0A: Checking samtools FASTA index"
echo "============================================================"

if [[ ! -s "$FAIDX" ]]; then
    echo "Creating FASTA index:"
    echo "$FAIDX"

    samtools faidx "$REF"
else
    echo "FASTA index already exists:"
    echo "$FAIDX"
fi

# ------------------------------------------------------------
# Step 0B: BWA index
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 0B: Checking BWA index"
echo "============================================================"

if [[ ! -s "${BWA_PREFIX}.bwt" ]]; then
    echo "Creating BWA index with prefix:"
    echo "$BWA_PREFIX"

    bwa index \
        -a bwtsw \
        -p "$BWA_PREFIX" \
        "$REF"
else
    echo "BWA index already exists:"
    echo "${BWA_PREFIX}.bwt"
fi

# ------------------------------------------------------------
# Step 1A: Map R1 separately
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 1A: Mapping Hi-C R1"
echo "============================================================"
echo "Start: $(date)"

bwa mem \
    -5SP \
    -t "$CPU" \
    "$BWA_PREFIX" \
    "$HIC_R1" |
    samtools view \
        -@ "$CPU" \
        -b \
        -o "${RAW_DIR}/${SAMPLE}_1.bam" \
        -

echo "R1 mapping finished: $(date)"

# ------------------------------------------------------------
# Step 1B: Map R2 separately
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 1B: Mapping Hi-C R2"
echo "============================================================"
echo "Start: $(date)"

bwa mem \
    -5SP \
    -t "$CPU" \
    "$BWA_PREFIX" \
    "$HIC_R2" |
    samtools view \
        -@ "$CPU" \
        -b \
        -o "${RAW_DIR}/${SAMPLE}_2.bam" \
        -

echo "R2 mapping finished: $(date)"

# ------------------------------------------------------------
# Step 2A: Filter R1 alignments at 5-prime end
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 2A: Filtering R1 at 5-prime end"
echo "============================================================"

samtools view \
    -@ "$CPU" \
    -h \
    "${RAW_DIR}/${SAMPLE}_1.bam" |
    perl "$FILTER" |
    samtools view \
        -@ "$CPU" \
        -b \
        -o "${FILT_DIR}/${SAMPLE}_1.bam" \
        -

# ------------------------------------------------------------
# Step 2B: Filter R2 alignments at 5-prime end
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 2B: Filtering R2 at 5-prime end"
echo "============================================================"

samtools view \
    -@ "$CPU" \
    -h \
    "${RAW_DIR}/${SAMPLE}_2.bam" |
    perl "$FILTER" |
    samtools view \
        -@ "$CPU" \
        -b \
        -o "${FILT_DIR}/${SAMPLE}_2.bam" \
        -

# ------------------------------------------------------------
# Step 3A: Combine read pairs, filter MAPQ and sort
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 3A: Pairing reads and applying MAPQ filter"
echo "MAPQ threshold: $MAPQ_FILTER"
echo "============================================================"

perl "$COMBINER" \
    "${FILT_DIR}/${SAMPLE}_1.bam" \
    "${FILT_DIR}/${SAMPLE}_2.bam" \
    samtools \
    "$MAPQ_FILTER" |
    samtools view \
        -b \
        -t "$FAIDX" \
        - |
    samtools sort \
        -@ "$CPU" \
        -T "${TMP_DIR}/${SAMPLE}.sort" \
        -o "${TMP_DIR}/${SAMPLE}.bam" \
        -

# ------------------------------------------------------------
# Step 3B: Add or replace read groups
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 3B: Adding read groups"
echo "============================================================"

java \
    -Xmx16G \
    -Djava.io.tmpdir="$TMP_DIR" \
    -jar "${EBROOTPICARD}/picard.jar" \
    AddOrReplaceReadGroups \
    INPUT="${TMP_DIR}/${SAMPLE}.bam" \
    OUTPUT="${PAIR_DIR}/${SAMPLE}.bam" \
    ID="$SAMPLE" \
    LB="$SAMPLE" \
    SM="$LABEL" \
    PL=ILLUMINA \
    PU="$SAMPLE" \
    SORT_ORDER=coordinate \
    VALIDATION_STRINGENCY=LENIENT

# ------------------------------------------------------------
# Step 4: Remove PCR/optical duplicates
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 4: Marking and removing duplicates"
echo "============================================================"

java \
    -Xmx48G \
    -XX:-UseGCOverheadLimit \
    -Djava.io.tmpdir="$TMP_DIR" \
    -jar "${EBROOTPICARD}/picard.jar" \
    MarkDuplicates \
    INPUT="${PAIR_DIR}/${SAMPLE}.bam" \
    OUTPUT="${REP_DIR}/${REP_LABEL}.bam" \
    METRICS_FILE="${REP_DIR}/metrics.${REP_LABEL}.txt" \
    TMP_DIR="$TMP_DIR" \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_DUPLICATES=true

# ------------------------------------------------------------
# Step 5: Index final BAM
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 5: Indexing final deduplicated BAM"
echo "============================================================"

samtools index \
    -@ "$CPU" \
    "${REP_DIR}/${REP_LABEL}.bam"

# ------------------------------------------------------------
# Step 6: Arima mapping statistics
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Step 6: Calculating mapping statistics"
echo "============================================================"

perl "$STATS" \
    "${REP_DIR}/${REP_LABEL}.bam" \
    > "${REP_DIR}/${REP_LABEL}.bam.stats"

# Additional samtools statistics
samtools flagstat \
    -@ "$CPU" \
    "${REP_DIR}/${REP_LABEL}.bam" \
    > "${REP_DIR}/${REP_LABEL}.bam.flagstat.txt"

samtools stats \
    -@ "$CPU" \
    "${REP_DIR}/${REP_LABEL}.bam" \
    > "${REP_DIR}/${REP_LABEL}.bam.samtools_stats.txt"

# ------------------------------------------------------------
# Final output information
# ------------------------------------------------------------

echo
echo "============================================================"
echo "Hi-C mapping pipeline completed successfully"
echo "============================================================"
echo "End time: $(date)"
echo
echo "Final deduplicated BAM:"
echo "${REP_DIR}/${REP_LABEL}.bam"
echo
echo "BAM index:"
echo "${REP_DIR}/${REP_LABEL}.bam.bai"
echo
echo "Arima statistics:"
echo "${REP_DIR}/${REP_LABEL}.bam.stats"
echo
echo "Picard duplicate metrics:"
echo "${REP_DIR}/metrics.${REP_LABEL}.txt"
echo
echo "samtools flagstat:"
echo "${REP_DIR}/${REP_LABEL}.bam.flagstat.txt"
echo "============================================================"
