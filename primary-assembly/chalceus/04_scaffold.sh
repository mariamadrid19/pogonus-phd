#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name scaffold
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o scaff.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

# Check completeness with minibusco
compleasm.py run -a purged.fa -o results_purged/ -l coleoptera -t 72

# sort the contigs by size
seqkit sort --by-length --reverse purged.fa > P_chalceus_sorted.fa

conda deactivate
conda activate ragtag

# obtain the first 20 scaffolds (20 largest ones) of the old reference genome
awk '/^>/ { if (count++ >= 20) exit } { print }'  old_ref.fa  > old_ref_20.fa 

# scaffold the contigs based on the old assembly (linkage groups)
ragtag.py scaffold old_ref.fa P_chalceus_sorted.fa -t 32 -o P_chalceus_scaffolds

# Check completeness with minibusco
compleasm.py run -a P_chalceus_scaffolds/ragtag.scaffold.fasta -o results_scaffolds/ -l coleoptera -t 72
