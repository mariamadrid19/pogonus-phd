#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=ont_filter
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o ont_filter.%j.out

seqtk seq -a GC157812_filtered.fastq > GC157812_filtered.fasta

module load Jellyfish/2.2.10-intel-2018a

#this will count kmers on the fasta file, takes a long time ( k = 17 is typical for long reads)
jellyfish count -m 17 -s 100M -t 32 -C GC157812_filtered.fasta

#this makes the k-mer histogram (can be seen on the .out file, DO NOT DELETE IT)
jellyfish histo -t 32 mer_counts.jf > mer_counts.histo

#this generates an incredibly big file! be careful with storage
jellyfish dump mer_counts.jf > mer_counts_dumps.fa

jellyfish info mer_counts.jf
