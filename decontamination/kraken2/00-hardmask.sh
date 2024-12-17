#!/bin/bash -l
#SBATCH --cluster=genius 
#SBATCH --job-name hardmask
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --output=hardmask.%j.out
#SBATCH -A lp_svbelleghem 

sed '/^>/!y/atcgn/NNNNN/' sorted_prim_dud.fasta.masked > sorted_prim_dud_hardmasked.fasta
grep -v '^>' sorted_prim_dud.fasta.masked | grep -o '[atcgn]' | wc -l # 807961487
grep -v '^>' sorted_prim_dud_hardmasked.fasta | grep -o '[atcgn]' | wc -l # 0
grep -v '^>' sorted_prim_dud_hardmasked.fasta | grep -o 'N' | wc -l # 807961487
