#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --partition=bigmem
#SBATCH --job-name=omnic_mapping
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=200G
#SBATCH --time=72:00:00
#SBATCH -o omnic_mapping.%j.out
#SBATCH -A lp_svbelleghem

##############################################
# Omni-C MAPPING PIPELINE #
##############################################

# Load required conda environment
conda activate omnic

# VARIABLES
CORES=36
REF='P_chalceus_REF1.fa'
GENOME_FILE='P_chalceus_REF1.fa.genome'
FAI_FILE="P_chalceus_REF1.fa.fai"
TMPDIR='tmpdir'
R1='TTCCAAGG-CCTTGTAG_R1.fastq.gz'
R2='TTCCAAGG-CCTTGTAG_R2.fastq.gz'
STATS='mapping_stats.txt'
PAIRS='mapped.pairs'
BAM='mapped.PT.bam'

# STEP 0: Index reference and prepare genome file
echo "### Step 0: Indexing reference"

# Check and create FASTA index (.fai)
if [ -f "$FAI_FILE" ]; then
    echo "FASTA index $FAI_FILE already exists. Skipping samtools faidx."
else
    echo "Creating FASTA index with samtools faidx..."
    samtools faidx "$REF"
fi

# Check and create genome file
if [ -f "$GENOME_FILE" ]; then
    echo "Genome file $GENOME_FILE already exists. Skipping cut."
else
    echo "Generating genome file..."
    cut -f1,2 "$FAI_FILE" > "$GENOME_FILE"
fi

# Check if BWA index exists (BWA creates several index files with various extensions)
if [ -f "${REF}.bwt" ]; then
    echo "BWA index files already exist. Skipping bwa index."
else
    echo "Creating BWA index..."
    bwa index "$REF"
fi

# Create temp directory if needed
[ -d $TMPDIR ] || mkdir -p $TMPDIR

# STEP 1: Mapping, parsing, sorting, deduplication, and conversion to BAM
echo "### Step 1: Mapping with bwa and processing with pairtools"

bwa mem -5SP -T0 -t $CORES $REF $R1 $R2 | \
pairtools parse \
  --min-mapq 40 \
  --walks-policy 5unique \
  --max-inter-align-gap 30 \
  --nproc-in $CORES \
  --nproc-out $CORES \
  --chroms-path $GENOME_FILE | \
pairtools sort \
  --tmpdir=$TMPDIR \
  --nproc $CORES | \
pairtools dedup \
  --nproc-in $CORES \
  --nproc-out $CORES \
  --mark-dups \
  --output-stats $STATS | \
pairtools split \
  --nproc-in $CORES \
  --nproc-out $CORES \
  --output-pairs $PAIRS \
  --output-sam - | \
samtools view -bS -@ $CORES - | \
samtools sort -@ $CORES -o $BAM

samtools index $BAM

# STEP 2: Mark duplicates in the BAM file using Picard
module load picard/2.18.23-Java-1.8.0_171

PICARD_OUT='mapped.marked.bam'
METRICS='mapped.dup_metrics.txt'

java -Xmx32G -Djava.io.tmpdir=$TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates \
  I=$BAM \
  O=$PICARD_OUT \
  M=$METRICS \
  TMP_DIR=$TMPDIR \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=LENIENT \
  REMOVE_DUPLICATES=true

samtools index $PICARD_OUT

echo "### Omni-C Mapping pipeline completed successfully"
