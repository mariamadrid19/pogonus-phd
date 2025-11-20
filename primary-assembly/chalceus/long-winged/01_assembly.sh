#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --partition=bigmem_long
#SBATCH --job-name Pog_hifiasm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=168:00:00
#SBATCH -o Pog_hifiasm.%j.out
#SBATCH -A lp_svbelleghem

#primary assembly
hifiasm -o Pchalceus_LW.asm --n-hap 2 --hom-cov 42 -t 24 --h1 GC170075_R1.fastq.gz --h2 GC170075_R2.fastq.gz raw_reads_LW_PacBio.fasta

#gfa_to_fasta
awk '/^S/{print ">"$2;print $3}' Pchalceus_LW.asm.hic.p_ctg.gfa > Pchalceus_LW.asm.hic.p_ctg.fa
