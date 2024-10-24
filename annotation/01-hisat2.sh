#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name hisat2
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH --output=hisat2.%j.out
#SBATCH -A lp_svbelleghem

conda activate isoseq

isoseq cluster POG_larveIsoSeq.demux.bam POG_unpolished_reads.demux.bam --verbose

samtools fasta POG_unpolished_reads.demux.bam POG_unpolished_reads.demux.fasta

module load HISAT2/2.1.0-intel-2018a
module load SAMtools/1.16.1-GCC-11.3.0

hisat2-build sorted_prim_dud.fasta sorted_prim_dud

hisat2 -f -x sorted_prim_dud -U POG_unpolished_reads.demux.fasta -S Dud_mapped_IsoSeq.sam

samtools view -h -q 30 -b Dud_mapped_IsoSeq.sam > Dud_mapped_IsoSeq_filtered.bam
samtools sort Dud_mapped_IsoSeq_filtered.bam > Dud_mapped_IsoSeq.filtered.bam
samtools index Dud_mapped_IsoSeq.filtered.bam

isoseq collapse --do-not-collapse-extra-5exons Dud_mapped_IsoSeq.filtered.bam collapsed_isoseq_dud_reads.gff
