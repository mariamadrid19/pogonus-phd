#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name Plit_hifiasm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o Plit_hifiasm.%j.out
#SBATCH -A lp_svbelleghem

hifiasm -o Pogonus_litt_l3.asm --primary --n-hap 2 --hom-cov 148 -l 3 -t 32 --h1 CTTGTCGA-GAACATCG_R1.fastq.gz --h2 CTTGTCGA-GAACATCG_R2.fastq.gz bc2042.fastq.gz

awk '$1 == "S" {print ">"$2"\n"$3}' Pogonus_litt_l3.asm.hic.p_ctg.gfa > Pogonus_litt_l3.asm.hic.p_ctg.fa

conda activate thesis 
assembly-stats Pogonus_litt_l3.asm.hic.p_ctg.fa
