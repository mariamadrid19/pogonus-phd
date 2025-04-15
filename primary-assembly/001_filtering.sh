#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=ont_filter
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=32G 
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o ont_filter.%j.out

# filter ONT long reads based on a minimum length of 5kb and the best 90% (based on quality) of those reads
filtlong --min_length 2000 --keep_percent 90 GC157812.fastq.gz | gzip >  GC157812_filtered.fastq

filtlong --min_length 2000 --keep_percent 90 GC157813.fastq.gz | gzip >  GC157813_filtered.fastq
