#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name softmask
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o softmask.%j.out
#SBATCH -A lp_svbelleghem


#to correctly install and run EDTA, had to install EDTA from yml file using mamba (https://github.com/oushujun/EDTA)
#and then had to download and make NINJA from source zip file (https://github.com/TravisWheelerLab/NINJA/releases/tag/0.98-cluster_only)
conda activate EDTA
 
BuildDatabase -name dudRM sorted_prim_dud.fasta
#need to tell RepeatModeler where the folder that holds NINJA is located
RepeatModeler -ninja_dir /data/leuven/357/vsc35707/NINJA-0.98-cluster_only/NINJA -database dudRM -pa 24 -LTRStruct >& run.out 
RepeatMasker -lib dudRM-families.fa -xsmall -pa 24 sorted_prim_dud.fasta
