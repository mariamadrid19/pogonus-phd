#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name BLAST
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o blast.%j.out
#SBATCH -A lp_svbelleghem

module load BLAST+/2.13.0-gompi-2022a

export PATH=${HOME}/edirect:${PATH}

sh get_species_taxids.sh -t 2 > bacterial.ids

sh get_species_taxids.sh -t 4751 > fungal.ids

sh get_species_taxids.sh -t 2157 > archaeal.ids

blastn -db nt \
       -query sorted_prim_dud.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -evalue 1e-5 \
       -num_threads 36 \
       -taxidlist bacterial.ids, fungal.ids, archaeal.ids \
       -out sorted_prim_dud.ncbi.blastn.out


