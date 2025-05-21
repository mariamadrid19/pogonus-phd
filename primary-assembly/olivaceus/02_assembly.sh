#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name Polv_hifiasm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=150:00:00
#SBATCH -o Polv_hifiasm.%j.out
#SBATCH -A lp_svbelleghem

# Assemble only 30% of the P. olivaceus reads
#hifiasm -o Pogonus_olivaceus.asm --n-hap 2 -t 36 P_olivaceus_HiFi_top30.fasta

#gfa_to_fasta
awk '/^S/{print ">"$2;print $3}' Pogonus_olivaceus.asm.p_ctg.gfa > Pogonus_olivaceus.asm.p_ctg.fa
