#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=ont_filter
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=32G 
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o ont_filter.%j.out

filtlong --min_length 1000 --min_mean_q 10 GC157812.fastq.gz | gzip >  GC157812_filtered.fastq

filtlong --min_length 1000 --min_mean_q 10 GC157813.fastq.gz | gzip >  GC157813_filtered.fastq
