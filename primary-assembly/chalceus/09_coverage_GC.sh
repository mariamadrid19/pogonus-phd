#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=cov_GC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=8:00:00
#SBATCH -o cov_GC.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis
module load BEDTools/2.27.1-intel-2018a

UNITIGS=/scratch/leuven/357/vsc35707/chalceus/03-T2T/Pogonus_T2T.asm.hic.p_utg.fa
INDEX=$UNITIGS.fai
RAWREADS=/scratch/leuven/357/vsc35707/chalceus/pacbio/GC157810.fasta
THREADS=16

# Map reads to unitigs and get sorted BAM
minimap2 -t $THREADS -x map-pb -a $UNITIGS $RAWREADS | samtools view -b - | samtools sort -@ 8 -o pacbio_mapped.bam

# Index the BAM
samtools index pacbio_mapped.bam

# Bedtools coverage per contig
bedtools genomecov -ibam pacbio_mapped.bam -g $INDEX -dz > coverage.txt

# Compute mean coverage per contig
awk '{cov[$1]+=$3; len[$1]++} END {for (i in cov) print i, cov[i]/len[i]}' coverage.txt > mean_coverage_per_contig.tsv

# GC content and length per contig
seqkit fx2tab -n -g $UNITIGS > gc_content.tsv

# Merge together
paste gc_content.tsv mean_coverage_per_contig.tsv > gc_vs_coverage.tsv
