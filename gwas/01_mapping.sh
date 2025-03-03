#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name bwa 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=48:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o map_RADTAGS.%j.out
#SBATCH --array=1-147

cd /scratch/leuven/357/vsc35707/GWAS/

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load BWA/0.7.17-foss-2018a
module load SAMtools/1.18-GCC-12.3.0
module load picard/2.18.23-Java-1.8.0_171
module load minimap2/2.26-GCCcore-12.3.0

echo "================="

# Sample IDs (146 samples)
samples=(Da_022 Da_029 Db_094 Db_098 Db_100 Db_101 Db_102 Db_103 GC3b_071 GC3b_072 GC3b_073 GC3b_074 GC3b_075 GC3b_076 GC3b_077 GC3b_078 GC3b_079 \
GC3b_080 GC3b_081 GC3b_082 GC3b_083 GC3b_084 GC3b_085 GC3b_086 GC3b_087 GC3b_088 GC3b_089 GC3b_090 GC3b_091 GC3b_092 GC3b_093 GC3b_098 GP3_055 \
GP3_059 GP3_067 GP3_068 GP3_078 GP3_088 GP3b_047 GP3b_048 GP3b_055 GP3b_056 GP3b_057 GP3b_062 GP3b_063 GP3b_064 GP3b_066 GP3b_067 GP3b_068 GP3b_069 \
GP3b_076 GP3b_077 GP3b_080 GP3b_082 GP3b_084 GP3b_085 Na_034 Nb_002 Nb_006 Nb_007 Nb_015 Nb_019 Nb_025 Nb_031 Nb_033 Nb_038 Nb_062 Nb_095 Pc2311 \
Pc2314 Pc2316 Pc2317 Pc2318 Pc2321 Pc2323 Pc2324 Pc2325 Pc2331 Pc2332 Pc2336 Pc2337 Pc2338 Pc2339 Pc2340 Pc277 Pc297 Pc298 Pc299 Pc300 Pc316 Pc317 \
Pc318 Pc786 Pc787 Pc788 Pc801 Pc805 Pc821 Pc822 Pc823 PcAVE_003 PcAVE_004 PcAVE_005 PcAVE_006 PcAVE_018 PcAVE_019 PcAVE_025 PcAVE_026 PcAVE_031 PcAVE_032 \
PcAVE_033 PcAVE_034 PcAVE_035 PcAVE_036 PcAVE_037 PcAVE_038 PcNP_033 PcNP_034 PcNP_035 PcNP_036 PcNP_037 PcNP_038 PcNP_040 PcNP_041 PcNP_043 PcNP_044 PcNP_045 \
PcNP_046 Pc_DZ_001 Pc_DZ_002 Pc_DZ_004 Pc_DZ_005 Pc_DZ_007 Pc_DZ_008 Pc_DZ_009 Pc_DZ_011 Pc_DZ_012 Pc_DZ_013 Pc_DZ_016 Pc_DZ_018 Pc_DZ_022 Pc_DZ_024 Pc_DZ_027 Pc_DZ_030)

echo "${samples[ID]}"

# Some folder and file paths to use later
REF=/scratch/leuven/357/vsc35707/GWAS/sorted_prim_dud.fasta
REFNAME=dudPrim
BWAout=/scratch/leuven/357/vsc35707/GWAS/bams
FILE1=/scratch/leuven/357/vsc35707/GWAS/fixed/$(echo "${samples[ID]}")_R1.fastq.gz
FILE2=/scratch/leuven/357/vsc35707/GWAS/fixed/$(echo "${samples[ID]}")_R2.fastq.gz

# Map reads using bwa mem
bwa mem -t 20 -M $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[ID]}").$REFNAME.bam

# Alternative, map with minimap2
minimap2 -ax sr -t 20 $REF $FILE1 $FILE2 | samtools view -bS - > $BWAout/$(echo "${samples[ID]}").$REFNAME.bam

# Filter using samtools
samtools view -f 0x02 -q 20 -b $BWAout/$(echo "${samples[ID]}").$REFNAME.bam > $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam

# Sort using samtools
samtools sort $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam -o $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam

# Remove PCR duplicates
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
INPUT=$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam \
OUTPUT=$BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam \
REMOVE_DUPLICATES=true \
METRICS_FILE=$BWAout/$(echo "${samples[ID]}").$REFNAME.dup_metrics.txt \
ASSUME_SORTED=true

# Remove intermediate files
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.bam
rm $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.bam 

samtools index $BWAout/$(echo "${samples[ID]}").$REFNAME.filtered.sorted.nd.bam
