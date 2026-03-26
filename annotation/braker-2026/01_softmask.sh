#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name softmask
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=72:00:00
#SBATCH -o softmask.%j.out
#SBATCH -A lp_svbelleghem

conda activate repeats
 
BuildDatabase -name PcSW Pchalceus_SW.sorted.fasta 
RepeatModeler -ninja_dir /data/leuven/357/vsc35707/NINJA-1.00-cluster_only/bin -database PcSW -threads 24 -LTRStruct >& run.out 
RepeatMasker -lib PcSW-families.fa -xsmall -pa 24 Pchalceus_SW.sorted.fasta 
