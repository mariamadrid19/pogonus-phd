#! /bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name map_contigs
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00
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
LABEL='Pogonus_gilvipes'
IN_DIR='/scratch/leuven/357/vsc35707/gilvipes/'
REF='/scratch/leuven/357/vsc35707/gilvipes/purged.fa'
FAIDX="${REF}.fai"
PREFIX='Pogonus_gilvipes'
RAW_DIR='/scratch/leuven/357/vsc35707/gilvipes/bams'
FILT_DIR='/scratch/leuven/357/vsc35707/gilvipes/filtered_bams'
FILTER='/scratch/leuven/357/vsc35707/gilvipes/filter_five_end.pl'
COMBINER='/scratch/leuven/357/vsc35707/gilvipes/two_read_bam_combiner.pl'
STATS='/scratch/leuven/357/vsc35707/gilvipes/get_stats.pl'
TMP_DIR='/scratch/leuven/357/vsc35707/gilvipes/temporary_files'
PAIR_DIR='/scratch/leuven/357/vsc35707/gilvipes/paired_bams'
REP_DIR='/scratch/leuven/357/vsc35707/gilvipes/deduplicated_files'
REP_LABEL=${LABEL}_r
MERGE_DIR='/scratch/leuven/357/vsc35707/gilvipes/final_merged_alignments'
MAPQ_FILTER=10
CPU=36

#load modules
module load picard/2.18.23-Java-1.8.0_171 #to run picard, use java -jar $EBROOTPICARD/picard.jar
conda activate thesis

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
bwa mem -5SP -t $CPU $REF ${SRA}_R1.fastq | samtools view -bS - > $RAW_DIR/${SRA}_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -5SP -t $CPU $REF ${SRA}_R2.fastq | samtools view -bS - > $RAW_DIR/${SRA}_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h $RAW_DIR/${SRA}_1.bam | perl $FILTER | samtools view -Sb - > $FILT_DIR/${SRA}_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h $RAW_DIR/${SRA}_2.bam | perl $FILTER | samtools view -Sb - > $FILT_DIR/${SRA}_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/${SRA}_1.bam $FILT_DIR/${SRA}_2.bam samtools $MAPQ_FILTER | samtools view -bS -t $FAIDX - | samtools sort -@ $CPU -o $TMP_DIR/$SRA.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups \
INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=PACBIO PU=none

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt \
TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"
