#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name eviann
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=12:00:00
#SBATCH -o eviann.%j.out
#SBATCH -A lp_svbelleghem

cd rnaseq/
paste <(ls $PWD/*_R1.fastq) <(ls $PWD/*_R2.fastq) > paired.txt

cd ../

eviann.sh -t 24 -g sorted_prim_dud.fasta -r paired.txt -p sequence.fasta -m 500000
