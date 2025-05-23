#! /bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name map_contigs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --time=48:00:00
#SBATCH -o map_contigs.%j.out
#SBATCH -A lp_svbelleghem

# This pipeline is modified from ArimaGenomics (arima_mapping_pipeline.sh)
# https://github.com/ArimaGenomics/mapping_pipeline/blob/master/arima_mapping_pipeline.sh

##############################################
# ARIMA GENOMICS MAPPING PIPELINE 07/26/2023 #
##############################################

# Below find the commands used to map HiC data to the contigs (output from hifiasm)
# This bash script will map one paired end HiC dataset (R1 & R2) to the contigs produced by hifiasm and purged by purge_dups (purged.fa)

##########################################
# Commands #
##########################################

SRA='CTTGTCGA-GAACATCG'
LABEL='Pogonus_littoralis'
IN_DIR='/scratch/leuven/357/vsc35707/littoralis/final_assembly'
REF='/scratch/leuven/357/vsc35707/littoralis/final_assembly/purged.fa'
FAIDX="${REF}.fai"
PREFIX='Pogonus_littoralis'
RAW_DIR='/scratch/leuven/357/vsc35707/littoralis/final_assembly/bams'
FILT_DIR='/scratch/leuven/357/vsc35707/littoralis/final_assembly/filtered_bams'
FILTER='/scratch/leuven/357/vsc35707/littoralis/final_assembly/filter_five_end.pl'
COMBINER='/scratch/leuven/357/vsc35707/littoralis/final_assembly/two_read_bam_combiner.pl'
STATS='/scratch/leuven/357/vsc35707/littoralis/final_assembly/get_stats.pl'
TMP_DIR='/scratch/leuven/357/vsc35707/littoralis/final_assembly/temporary_files'
PAIR_DIR='/scratch/leuven/357/vsc35707/littoralis/final_assembly/paired_bams'
REP_DIR='/scratch/leuven/357/vsc35707/littoralis/final_assembly/deduplicated_files'
REP_LABEL=${LABEL}_r
MERGE_DIR='/scratch/leuven/357/vsc35707/littoralis/final_assembly/final_merged_alignments'
MAPQ_FILTER=10
CPU=36

#load modules
module load picard/2.18.23-Java-1.8.0_171 #to run picard, use java -jar $EBROOTPICARD/picard.jar
module load SAMtools/0.1.20-GCC-12.3.0
module load BWA/0.7.17-GCC-10.3.0
mamba activate pairtools

echo "### Step 0: Check output directories' existence & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
if [ ! -e "${PREFIX}.bwt" ]; then
    bwa index -a bwtsw -p $PREFIX $REF
fi

echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -5SP -t $CPU $REF ${SRA}_R1.fastq.gz | samtools view -bS - > $RAW_DIR/${SRA}_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -5SP -t $CPU $REF ${SRA}_R2.fastq.gz | samtools view -bS - > $RAW_DIR/${SRA}_2.bam

echo "### Step 2.A: Sort BAMs by read name (1st)"
samtools sort -n -@ $CPU -o $RAW_DIR/${SRA}_1.qsort.bam $RAW_DIR/${SRA}_1.bam

echo "### Step 2.A: Sort BAMs by read name (2nd)"
samtools sort -n -@ $CPU -o $RAW_DIR/${SRA}_2.qsort.bam $RAW_DIR/${SRA}_2.bam

echo "### Step 3: Merge mates and extract pairs"
pairtools parse \
    --no-flip \
    --min-mapq 0 \
    --nproc-in $CPU \
    --nproc-out $CPU \
    --drop-sam \
    --output $PAIR_DIR/${SRA}.pairs \
    --output-stats $PAIR_DIR/${SRA}.parse_stats.txt \
    $RAW_DIR/${SRA}_1.qsort.bam $RAW_DIR/${SRA}_2.qsort.bam

echo "### Step 4: Sort and deduplicate pairs"
pairtools sort -p $PAIR_DIR/${SRA}.pairs --output $PAIR_DIR/${SRA}.sorted.pairs
pairtools dedup --output-stats $PAIR_DIR/${SRA}.dedup_stats.txt $PAIR_DIR/${SRA}.sorted.pairs > $PAIR_DIR/${SRA}.dedup.pairs

echo "### Step 5: Filter for species-specific high-quality pairs"
awk -F'\t' '{
    if ($1 ~ /^#/) {print; next}
    if ($8 == "UU" && $12 >= 20 && $13 >= 20 &&
        $0 ~ /NM:i:([0-2])\t.*\t.*\t.*\t.*\t.*NM:i:([0-2])/) {
        print
    }
}' $PAIR_DIR/${SRA}.dedup.pairs > $FILT_DIR/${SRA}.filtered.pairs


# Convert pairs to BAM
pairtools split --output-pairs $FILT_DIR/${SRA}.split.pairs \
                --output-pairs-sam $FILT_DIR/${SRA}.filtered.sam \
                $FILT_DIR/${SRA}.filtered.pairs

# Convert to BAM, sort, and index
samtools view -Sb $FILT_DIR/${SRA}.filtered.sam | \
samtools sort -@ $CPU -o $PAIR_DIR/${SRA}.bam
samtools index $PAIR_DIR/${SRA}.bam

# Step 6: Add read group
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
    INPUT=$PAIR_DIR/${SRA}.bam \
    OUTPUT=$PAIR_DIR/${SRA}.rg.bam \
    ID=$SRA LB=$SRA SM=$LABEL PL=PACBIO PU=none

# Step 7: Mark duplicates
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    INPUT=$PAIR_DIR/${SRA}.rg.bam \
    OUTPUT=$REP_DIR/$REP_LABEL.bam \
    METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt \
    TMP_DIR=$TMP_DIR \
    ASSUME_SORTED=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    REMOVE_DUPLICATES=TRUE

samtools index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"
