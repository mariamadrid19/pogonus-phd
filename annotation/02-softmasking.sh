#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name softmask
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o softmask.%j.out
#SBATCH -A lp_svbelleghem

conda activate EDTA2
 
BuildDatabase -name dudRM sorted_prim_dud.fasta
RepeatModeler -database dudRM -threads 24 -LTRStruct >& run.out
RepeatMasker -lib dudRM-families.fa -xsmall -pa 24 sorted_prim_dud.fasta
