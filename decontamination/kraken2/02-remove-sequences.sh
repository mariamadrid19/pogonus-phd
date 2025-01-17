#!/bin/bash -l
#SBATCH --cluster=genius 
#SBATCH --job-name kraken2
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --output=kraken2.%j.out
#SBATCH -A lp_svbelleghem 

#Cut columns: c1,c2,c3

#Filter: Cut on output of Cut with following condition: c1!='U'

#Cut columns: c2 from the output of Filtering

#Rename the output "sequences_to_remove.txt", it's sequence IDs to be removed 

#Remove sequences with gfastats

gfastats sorted_prim_dud.fasta --remove sequences_to_remove.txt -o dudzele_decont.fasta
