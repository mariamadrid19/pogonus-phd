#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name jf_reseq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -o jf_reseq.%j.out
#SBATCH -A lp_svbelleghem

module load Jellyfish/2.2.10-intel-2018a

# Estimate genome size based on short reads (Ilummina WGS, reseq), zipped
# Count k-mers with Jellyfish
zcat Bar2_01_R1.fq.gz Bar2_01_R2.fq.gz | jellyfish count -C -m 21 -s 10G -t 24 --min-quality-char=? -o mer_counts_reseq.jf
# This .histo file can then be used as input for GenomeScope2
jellyfish histo -t 32 mer_counts_reseq.jf > mer_counts_reseq.histo
cp mer_counts_reseq.histo $VSC_DATA
