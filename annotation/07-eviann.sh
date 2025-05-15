#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name eviann
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH -o eviann.%j.out
#SBATCH -A lp_svbelleghem

mamba activate entrez
# use esearch and efetch to download all the proteins for insects available in NCBI
esearch -db protein -query "txid50557[Organism:exp]" | efetch -format fasta > insecta_proteins.fasta

# generate a txt file with all the names of the RNAseq reads needed (they should all be in a directory called rnaseq)
cd rnaseq/
paste <(ls $PWD/*_R1.fastq) <(ls $PWD/*_R2.fastq) > paired_mixed.txt

cd ../

# run EviAnn with the proteins and the RNA sequences
eviann.sh -t 24 -g sorted_prim_dud.fasta -r paired_mixed.txt -p insecta_proteins.fasta -m 500000 --lncrnamintpm 1
