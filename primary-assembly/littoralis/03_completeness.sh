#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name compleasm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH -o compleasm.%j.out
#SBATCH -A lp_svbelleghem

conda deactivate
compleasm.py run -a Pogonus_litt_l3.asm.hic.p_ctg.fa -o results_l3/ -l coleoptera -t 32
