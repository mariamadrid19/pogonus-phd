#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name scaffold
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o scaff.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

#Check completeness with minibusco
compleasm.py run -a purged.fa -o results_purged/ -l coleoptera -t 72

#sort the contigs by size
seqkit sort --by-length --reverse purged.fa > P_littoralis_sorted.fa

#obtain the first 50 scaffolds (50 largest ones)
awk '/^>/ { if (count++ >= 50) exit } { print }' P_littoralis_sorted.fa > P_littoralis_sorted_50.fa

conda deactivate
conda activate ragtag

# scaffold the contigs based on the old assembly (linkage groups)
ragtag.py scaffold -f -w 1000 -q 60 old_ref.fa P_littoralis_sorted_50.fa

compleasm.py run -a P_littoralis_sorted_50.fa -o results_scaffolds/ -l coleoptera -t 72
