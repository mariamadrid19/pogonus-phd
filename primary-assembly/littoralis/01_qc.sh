#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name Plit_jf
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o Plit_jf.%j.out
#SBATCH -A lp_svbelleghem

module load Jellyfish/2.2.10-intel-2018a

#this will count kmers on the fasta file, takes a long time
jellyfish count -m 31 -s 100M -t 32 -C GC157811.fasta

#this makes the k-mer histogram (can be seen on the .out file, DO NOT DELETE IT)
jellyfish histo -t 32 mer_counts.jf > mer_counts.histo

mv mer_counts.histo Plit_k31.histo
rm mer_counts.jf

cp Plit_k31.histo $VSC_DATA # use this hisogram to generate the GenomeScope estimate
