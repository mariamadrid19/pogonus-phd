#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name compleasm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -o compleasm.%j.out
#SBATCH -A lp_svbelleghem

#conda install pandas
#conda deactivate
#run outside of any conda env (won't run)

# Run compleasm
compleasm_kit/compleasm.py download archaea
compleasm_kit/compleasm.py download bacteria
compleasm_kit/compleasm.py download fungi

compleasm_kit/compleasm.py run -a sorted_prim_dud.fasta -o archaea/ -l archaea -t 24
compleasm_kit/compleasm.py run -a sorted_prim_dud.fasta -o fungi/ -l fungi -t 24
compleasm_kit/compleasm.py run -a sorted_prim_dud.fasta -o bacteria/ -l bacteria -t 24
