#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name isoseq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH -o isoseq.%j.out
#SBATCH -A lp_svbelleghem

conda activate isoseq

isoseq collapse --do-not-collapse-extra-5exons Dud_mapped_IsoSeq.filtered.bam collapsed_isoseq.gff
