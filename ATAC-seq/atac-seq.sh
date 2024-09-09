#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name ATAC_test 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8 
#SBATCH --time=24:10:00 
#SBATCH -A lp_svbelleghem
 
module load Bowtie2/2.3.4.3-foss-2018a
module load SAMtools/1.16.1-GCC-11.3.0
module load Java/11.0.20
module load BEDTools/2.27.1-GCCcore-6.4.0
 
ID=$((SLURM_ARRAY_TASK_ID -1))
 
samples=(GC154460)
 
module load Trimmomatic/0.39-Java-1.8.0_192
 
cd /scratch/leuven/357/vsc35707/Pogonus_ATACseq 

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
$(echo "${samples[ID]}")_R1.fastq.gz $(echo "${samples[ID]}")_R2.fastq.gz \
$(echo "${samples[ID]}")_R1_trimmed.fastq.gz $(echo "${samples[ID]}")_R1_unpaired.fastq.gz \
$(echo "${samples[ID]}")_R2_trimmed.fastq.gz $(echo "${samples[ID]}")_R2_unpaired.fastq.gz \
ILLUMINACLIP:/scratch/leuven/357/vsc35707/Pogonus_ATACseq/TruSeq3-PE.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:50 HEADCROP:10
 
cd /scratch/leuven/357/vsc35707/Pogonus_ATACseq
 
REF=/scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta
 
samtools faidx /scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta; cut -f1,2 /scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta.fai > /scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta.sizes
 
SIZES=/scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta.sizes
 
bowtie2-build /scratch/leuven/357/vsc35707/Pogonus_ATACseq/sorted_prim_dud.fasta sorted_prim_dud
 
bowtie2 -t -k 2 -p 8 --local -x $REF \
-1 /scratch/leuven/357/vsc35707/Pogonus_ATACseq/$(echo "${samples[ID]}")_R1_trimmed.fastq.gz \
-2 /scratch/leuven/357/vsc35707/Pogonus_ATACseq/$(echo "${samples[ID]}")_R2_trimmed.fastq.gz \
samtools view -bS - > $(echo "${samples[ID]}").bam
 
samtools view -f 0x02 -q 20 -b $(echo "${samples[ID]}").bam > $(echo "${samples[ID]}").filtered.bam
 
samtools sort $(echo "${samples[ID]}").filtered.bam -o $(echo "${samples[ID]}").filtered.sorted.bam

module load picard/2.18.23-Java-1.8.0_171

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=$(echo "${samples[ID]}").filtered.sorted.bam \
OUTPUT=$(echo "${samples[ID]}").filtered.sorted.nd.bam \
METRICS_FILE=$(echo "${samples[ID]}")_dup_metrics.txt \
ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE
 
##
 
bedtools genomecov \
-ibam $(echo "${samples[ID]}").filtered.sorted.nd.bam -bg \
> $(echo "${samples[ID]}").filtered.sorted.nd.bdg
 
LC_COLLATE=C sort -k1,1 -k2,2n $(echo "${samples[ID]}").filtered.sorted.nd.bdg \
> $(echo "${samples[ID]}").filtered.sorted.nd.collate.bdg
 
conda activate bigwig
 
bedGraphToBigWig \
$(echo "${samples[ID]}").filtered.sorted.nd.collate.bdg \
$SIZES \
$(echo "${samples[ID]}").filtered.sorted.nd.collate.bw
 
