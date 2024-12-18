#!/bin/bash -l
#SBATCH --cluster=genius 
#SBATCH --job-name kraken2
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --output=kraken2.%j.out
#SBATCH -A lp_svbelleghem 

DBNAME=pogonus_dud

kraken2-build --standard --threads 32 --db $DBNAME

kraken2 --db $DBNAME --threads 32 sorted_prim_dud_hardmasked.fasta
