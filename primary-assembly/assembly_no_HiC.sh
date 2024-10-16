#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name nieu_no_hic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o assembly.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/nieuwpoort_pogonus/

#assemble the contigs using hifiasm, without any HiC data
hifiasm -o prim_nieu_no_hic.asm --n-hap 2 --hom-cov 30 -t 32 POG_Nieuwpoort_HiFi_reads.fasta

#gfa to fasta
awk '/^S/{print ">"$2;print $3}' prim_nieu_no_hic.asm.bp.p_ctg.gfa > prim_nieu_no_hic.asm.bp.p_ctg.fa

module load SAMtools/1.13-GCC-10.3.0
module load BWA/0.7.17-GCC-10.3.0

#indexing
samtools faidx prim_nieu_no_hic.asm.bp.p_ctg.fa && cut -f1,2 prim_nieu_no_hic.asm.bp.p_ctg.fa.fai > prim_nieu_no_hic.asm.bp.p_ctg.fa.genome && bwa index prim_nieu_no_hic.asm.bp.p_ctg.fa

#variables for the HiC mapping pipeline
SRA='GC143248_ACTCTCGA-TGGTACAG_S65'
LABEL='Pogonus_chalceus'
IN_DIR='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus'
REF='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/prim_nieu_no_hic.asm.bp.p_ctg.fa'
FAIDX='$REF.fai'
PREFIX='Pogonus_hifiasm.asm.hic.p_ctg'
RAW_DIR='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/bams'
FILT_DIR='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/filtered_bams'
FILTER='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/filter_five_end.pl'
COMBINER='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/two_read_bam_combiner.pl'
STATS='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/get_stats.pl'
TMP_DIR='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/temporary_files'
PAIR_DIR='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/paired_bams'
REP_DIR='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/deduplicated_files'
REP_LABEL=${LABEL}_r
MERGE_DIR='/scratch/leuven/357/vsc35707/nieuwpoort_pogonus/map_contigs/final_merged_alignments'
MAPQ_FILTER=10
CPU=12

conda activate thesis 

#load picard module
module load picard/2.18.23-Java-1.8.0_171

echo "### Step 0: Check output directories' existence & create them as needed"
[ -d $RAW_DIR ] || mkdir -p $RAW_DIR
[ -d $FILT_DIR ] || mkdir -p $FILT_DIR
[ -d $TMP_DIR ] || mkdir -p $TMP_DIR
[ -d $PAIR_DIR ] || mkdir -p $PAIR_DIR
[ -d $REP_DIR ] || mkdir -p $REP_DIR
[ -d $MERGE_DIR ] || mkdir -p $MERGE_DIR

echo "### Step 0: Index reference" # Run only once! Skip this step if you have already generated BWA index files
bwa index -a bwtsw -p $PREFIX $REF

echo "### Step 1.A: FASTQ to BAM (1st)"
bwa mem -t $CPU $REF $IN_DIR/${SRA}_R1.fastq | samtools view -@ $CPU -Sb - > $RAW_DIR/${SRA}_1.bam

echo "### Step 1.B: FASTQ to BAM (2nd)"
bwa mem -t $CPU $REF $IN_DIR/${SRA}_R2.fastq | samtools view -@ $CPU -Sb - > $RAW_DIR/${SRA}_2.bam

echo "### Step 2.A: Filter 5' end (1st)"
samtools view -h $RAW_DIR/${SRA}_1.bam | perl $FILTER | samtools view -Sb - > $FILT_DIR/${SRA}_1.bam

echo "### Step 2.B: Filter 5' end (2nd)"
samtools view -h $RAW_DIR/${SRA}_2.bam | perl $FILTER | samtools view -Sb - > $FILT_DIR/${SRA}_2.bam

echo "### Step 3A: Pair reads & mapping quality filter"
perl $COMBINER $FILT_DIR/${SRA}_1.bam $FILT_DIR/${SRA}_2.bam samtools $MAPQ_FILTER | samtools view -bS -t $FAIDX - | samtools sort -@ $CPU -o $TMP_DIR/$SRA.bam -

echo "### Step 3.B: Add read group"
java -Xmx4G -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups INPUT=$TMP_DIR/$SRA.bam OUTPUT=$PAIR_DIR/$SRA.bam ID=$SRA LB=$SRA SM=$LABEL PL=PACBIO PU=none

echo "### Step 4: Mark duplicates"
java -Xmx30G -XX:-UseGCOverheadLimit -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$PAIR_DIR/$SRA.bam OUTPUT=$REP_DIR/$REP_LABEL.bam METRICS_FILE=$REP_DIR/metrics.$REP_LABEL.txt TMP_DIR=$TMP_DIR ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

samtools index $REP_DIR/$REP_LABEL.bam

perl $STATS $REP_DIR/$REP_LABEL.bam > $REP_DIR/$REP_LABEL.bam.stats

echo "Finished Mapping Pipeline through Duplicate Removal"



cd /scratch/leuven/357/vsc35707/nieuwpoort_pogonus/

#scaffold the contigs using the primary assembly and the bam files (HiC reads mapped to the contigs)
yahs prim_nieu_no_hic.asm.bp.p_ctg.fa /map_contigs/deduplicated_files/Pogonus_chalceus_r.bam -q 100000 -l 30 -r 50000 -o nieu_scaffolds
