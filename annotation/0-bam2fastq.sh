#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name b2fq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH -o b2fq.%j.out
#SBATCH -A lp_svbelleghem

conda activate pbtk

bam2fastq -o out POG_IsoSeq_HiFi_demux POG_larveIsoSeq.demux.bam

cat POG_IsoSeq_HiFi_demux.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > POG_IsoSeq_HiFi_demux.fasta
