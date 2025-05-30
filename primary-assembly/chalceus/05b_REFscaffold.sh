#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name scaffold
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o scaff.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis
seqkit sort --by-length --reverse P_chalceus_REF1.fa  > P_chalceus_REF1_sorted.fa
awk '/^>/ { if (count++ >= 11) exit } { print }' old_ref.fa > 11_old_ref.fa

conda deactivate
conda activate ragtag

# scaffold the contigs based on the old assembly (linkage groups)
ragtag.py scaffold -u 11_old_ref.fa P_chalceus_REF1_sorted.fa -t 32 -o REF1_scaffolds

seqkit sort --by-length --reverse REF1_scaffolds/ragtag.scaffold.fasta > REF1_scaffolds/ragtag.scaffold.sorted.fasta

samtools faidx REF1_scaffolds/ragtag.scaffold.fasta
samtools faidx REF1_scaffolds/ragtag.scaffold.sorted.fasta

# Count how many of the scaffolds were placed as scaffolds based on the old reference genome (RagTag scaffolds)
awk '$1 ~ /^CM00/ {sum += $2} END {print sum}' REF1_scaffolds/ragtag.scaffold.fasta.fai
awk '$1 ~ /_RagTag$/ {sum += $2} END {print sum}' REF1_scaffolds/ragtag.scaffold.fasta.fai

# Check completeness with minibusco
compleasm.py run -a REF1_scaffolds/ragtag.scaffold.fasta -o results_scaffolds/ -l coleoptera -t 72
