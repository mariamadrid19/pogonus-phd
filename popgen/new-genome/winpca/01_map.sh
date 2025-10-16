#!/bin/bash -l 
#SBATCH --cluster=genius
#SBATCH --job-name map
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=36
#SBATCH --time=60:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o map_all.%j.out
#SBATCH --array=1-103

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

conda activate variant_tools
# samtools, bwa, vcftools are here 

echo "================="

# Sample IDs (all samples, 103)
samples=(\
GC129388 GC129389 GC129390 GC129391 GC129392 \
GC129393 GC129394 GC129395 GC129396 GC129397 \
GC129398 GC129399 GC129400 GC129401 GC129402 \
GC129403 GC129404 GC129405 GC129406 GC129407 \
GC129408 GC129409 GC129410 GC129411 GC129412 \
GC129413 GC129414 GC129415 GC129416 GC129417 \
GC129418 GC129419 GC129420 GC129421 GC129422 \
GC129423 GC129424 GC129425 GC129426 GC129427 \
GC129428 GC129429 GC129430 GC129431 GC129432 \
GC129433 GC129434 GC129435 GC129437 GC129438 \
GC129439 GC129440 GC136078 GC136079 GC136080 \
GC136081 GC136082 GC136083 GC136084 GC136085 \
GC136086 GC136087 GC136088 GC136089 GC136090 \
GC136091 GC136092 GC136093 GC136094 GC136095 \
GC136096 GC136097 GC136098 GC136099 GC136100 \
GC136101 GC136102 GC136103 GC136104 GC136105 \
GC136106 GC136107 GC136108 GC136109 GC136110 \
GC136111 GC136112 GC136113 GC136114 GC136115 \
GC136116 GC136117 GC136118 GC136119 GC136120 \
GC136121 GC136122 GC136123 GC136124 GC136125 \
GC136126 GC136127 GC136128)

echo "${samples[ID]}"

# Some folder and file paths to use later
REF=/scratch/leuven/357/vsc35707/winpca/Pchalceus_SW.sorted.fasta
REFNAME=Pchal_Bar_SW
BWAout=/scratch/leuven/357/vsc35707/winpca/bams
FILE1=/scratch/leuven/357/vsc35707/winpca/reads/$(echo "${samples[ID]}")_R1.fastq.gz
FILE2=/scratch/leuven/357/vsc35707/winpca/reads/$(echo "${samples[ID]}")_R2.fastq.gz

mkdir -p $BWAout

echo "Mapping reads: ${samples[ID]}"
# Map reads using bwa mem
bwa mem -t 36 -M $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[ID]}").$REFNAME.bam

echo "Filtering reads: ${samples[ID]}"
# Filter using samtools
samtools view -f 0x02 -q 20 -b $BWAout/$(echo "${samples[ID]}").$REFNAME.bam > $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam

# Remove .bam if .filtered.bam is larger than 500 KB
if [ $(stat -c %s "$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam") -gt 512000 ]; then
    rm "$BWAout/$(echo "${samples[ID]}").$REFNAME.bam"
fi

echo "Sorting reads: ${samples[ID]}"
# Sort using samtools
samtools sort $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam -o $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

# Remove .filtered.bam if .filtered.sorted.bam is larger than 500 KB
if [ $(stat -c %s "$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam") -gt 512000 ]; then
    rm "$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam"
fi

module load picard/2.18.23-Java-1.8.0_171

echo "Removing PCR duplicates: ${samples[ID]}"
# Remove PCR duplicates
java -Djava.io.tmpdir=temp/ -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam O=$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.dedup.bam REMOVE_DUPLICATES=true M=$BWAout/$(echo "${samples[ID]}").$REFNAME.dup_metrics.txt ASSUME_SORTED=true

# Remove .filtered.sorted.bam if .dedup.bam is larger than 500 KB
if [ $(stat -c %s "$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.dedup.bam") -gt 512000 ]; then
    rm "$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam"
fi

samtools index $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.dedup.bam

echo "Finished with ${samples[ID]}"
