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
bam2fastq -o bc2042 hifi_reads/m84247_241127_161129_s3.hifi_reads.bc2042.bam

module load FastQC/0.12.1-Java-11

fastqc bc2041.fastq.gz -t 32 #SW
fastqc bc2042.fastq.gz -t 32 #LW 
