#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name qc
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=48:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o bams_qc.%j.out

cd /scratch/leuven/357/vsc35707/GWAS/bams

for bam in *.bam; do 
    echo "$bam"; 
    samtools flagstat "$bam"; 
done > mapping_stats.txt

for bam in *.bam; do 
    echo -n "$(basename "$bam" .bam): "; 
    samtools depth -a "$bam" | awk '{sum+=$3} END {if (NR>0) print sum/NR; else print "No coverage"}'; 
done > average_depths.txt
