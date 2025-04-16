#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name jellyfish
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH -o jellyfish.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis 

#this will convert the fastq.gz reads into fasta reads
seqtk seq -a bc2041.fastq.gz > GC157810.fasta # SW
seqtk seq -a bc2042.fastq.gz > GC157811.fasta # LW

module purge
module load Jellyfish/2.2.10-intel-2018a

#this will count kmers on the fasta file, takes a long time
jellyfish count -m 21 -s 100M -t 32 -C GC157810.fasta

#this makes the k-mer histogram (can be seen on the .out file, DO NOT DELETE IT)
jellyfish histo -t 32 mer_counts.jf > mer_counts.histo

#this generates an incredibly big file! be careful with storage
jellyfish dump mer_counts.jf > mer_counts_dumps.fa

jellyfish info mer_counts.jf
