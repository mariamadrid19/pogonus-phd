#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name ATAC_test 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:10:00 
#SBATCH -A lp_svbelleghem
 
module load Bowtie2/2.3.4.3-foss-2018a
module load SAMtools/1.16.1-GCC-11.3.0
module load Java/11.0.20
module load BEDTools/2.27.1-GCCcore-6.4.0
module load Trimmomatic/0.39-Java-1.8.0_192
module load picard/2.18.23-Java-1.8.0_171
 
cd /scratch/leuven/357/vsc35707/Pogonus_ATACseq 

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
GC154460_R1.fastq.gz GC154460_R2.fastq.gz \
GC154460_R1_trimmed.fastq.gz GC154460_R1_unpaired.fastq.gz \
GC154460_R2_trimmed.fastq.gz GC154460_R2_unpaired.fastq.gz \
ILLUMINACLIP:/scratch/leuven/357/vsc35707/Pogonus_ATACseq/TruSeq3-PE.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50 HEADCROP:10
  
REF=/scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta
 
samtools faidx sorted_prim_dud.fasta; cut -f1,2 sorted_prim_dud.fasta.fai > sorted_prim_dud.fasta.sizes
 
SIZES=/scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta.sizes
 
bowtie2-build /scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta sorted_prim_dud
 
bowtie2 -t -k 2 -p 8 --local -x sorted_prim_dud.fasta -1 GC154460_R1_trimmed.fastq.gz -2 GC154460_R2_trimmed.fastq.gz | samtools view -bS - > GC154460.bam
 
samtools view -f 0x02 -q 20 -b GC154460.bam > GC154460.filtered.bam
 
samtools sort GC154460.filtered.bam -o GC154460.filtered.sorted.bam

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=GC154460.filtered.sorted.bam \
OUTPUT=GC154460.filtered.sorted.nd.bam \
METRICS_FILE=GC154460_dup_metrics.txt \
ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE
 
bedtools genomecov -ibam GC154460.filtered.sorted.nd.bam -bg > GC154460.filtered.sorted.nd.bdg
 
LC_COLLATE=C sort -k1,1 -k2,2n GC154460.filtered.sorted.nd.bdg > GC154460.filtered.sorted.nd.collate.bdg
 
conda activate bigwig
 
bedGraphToBigWig GC154460.filtered.sorted.nd.collate.bdg $SIZES GC154460.filtered.sorted.nd.collate.bw
