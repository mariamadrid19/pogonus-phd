#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name sort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o sort.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis
seqkit sort --by-length --reverse purged.fa > P_littoralis_sorted.fa

conda deactivate
conda activate ragtag

# scaffold the contigs based on the old assembly (linkage groups)
ragtag.py scaffold old_ref.fa P_littoralis_sorted.fa -t 32 -o P_littoralis_scaffolds

compleasm.py run -a P_littoralis_scaffolds/ragtag.scaffold.fasta -o results_scaffolds/ -l coleoptera -t 72
