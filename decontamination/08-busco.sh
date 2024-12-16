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
compleasm.py download archaea
compleasm.py download bacteria
compleasm.py download fungi

compleasm.py run -a sorted_prim_dud.fasta -o archaea/ -l archaea -t 24
compleasm.py run -a sorted_prim_dud.fasta -o fungi/ -l fungi -t 24
compleasm.py run -a sorted_prim_dud.fasta -o bacteria/ -l bacteria -t 24


busco -i sorted_prim_dud.fasta -m genome -l archaea -c 24 -o archaea_busco/
busco -i sorted_prim_dud.fasta -m genome -l fungi -c 24 -o fungi_busco/
busco -i sorted_prim_dud.fasta -m genome -l bacteria -c 24 -o bacteria_busco/
