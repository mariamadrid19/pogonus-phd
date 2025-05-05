#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name pb_qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=10:00:00
#SBATCH -o pb_qc.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/chalceus/pacbio

conda activate thesis
bam2fastq -o bc2041 hifi_reads/m84247_241127_161129_s3.hifi_reads.bc2041.bam

module load FastQC/0.12.1-Java-11
fastqc bc2041.fastq.gz -t 32

# fastq to fasta
seqtk seq -a bc2041.fastq.gz > GC157810.fasta

module purge
module load Jellyfish/2.2.10-intel-2018a

#this will count kmers on the fasta file, takes a long time (k = 31 for repeatitive genomes)
jellyfish count -m 31 -s 100M -t 32 -C GC157810.fasta

#this makes the k-mer histogram, open it with GenomeScope2
jellyfish histo -t 32 mer_counts.jf > mer_counts.histo

module purge
module load NanoPlot/1.42.0-foss-2022a

# quality check of the reads
NanoPlot --fastq bc2041.fastq.gz -o P_chalceus_Nanoplot
