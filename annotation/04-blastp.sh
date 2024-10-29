#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name blastP
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o blastP.%j.out
#SBATCH -A lp_svbelleghem

module load BLAST+/2.13.0-gompi-2022a

#all the NCBI BLAST databases: https://ftp.ncbi.nlm.nih.gov/blast/db/

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/swissprot.tar.gz.md5

tar -xvzf swissprot.tar.gz
  
blastp -db swissprot -query ../augustus.hints.h.aa -outfmt 5 -max_target_seqs 10 -max_hsps 1 -out blastout.swiss.h.xml
