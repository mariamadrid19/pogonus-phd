#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name pb_qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH -o pb_qc.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis
bam2fastq -o raw_reads_LW_PacBio m84288_251113_172633_s2.hifi_reads.bc2001.bam

module load FastQC/0.12.1-Java-11
fastqc raw_reads_LW_PacBio.fastq.gz -t 32

# fastq to fasta
seqtk seq -a raw_reads_LW_PacBio.fastq.gz > raw_reads_LW_PacBio.fasta

#this will count kmers on the fasta file, takes a long time (k = 31 for highly repeatitive genomes)
jellyfish count -m 31 -s 100M -t 32 -C raw_reads_LW_PacBio.fasta

#this makes the k-mer histogram, open it with GenomeScope2
jellyfish histo -t 32 mer_counts.jf > mer_counts.histo

module purge
module load NanoPlot/1.42.0-foss-2022a

# quality check of the reads
NanoPlot --fastq raw_reads_LW_PacBio.fastq.gz -o P_chalceus_LW_Nanoplot
