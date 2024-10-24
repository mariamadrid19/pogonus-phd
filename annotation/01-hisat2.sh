#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name hisat2
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH --output=hisat2.%j.out
#SBATCH -A lp_svbelleghem

#Cluster IsoSeq reads, map these transcripts to the assembled genome, and collapse these mapped reads

conda activate isoseq

isoseq cluster POG_larveIsoSeq.demux.bam POG_unpolished_reads.demux.bam --verbose

gunzip POG_unpolished_reads.demux.hq.fasta.gz

module load HISAT2/2.1.0-intel-2018a

hisat2-build sorted_prim_dud.fasta sorted_prim_dud

hisat2 --dta -f -x sorted_prim_dud -U POG_unpolished_reads.demux.hq.fasta -S POG_mapped_RNA_dud.sam

module load SAMtools/1.16.1-GCC-11.3.0

samtools view -h -q 30 -b POG_mapped_RNA_dud.sam > POG_mapped_RNA_dud.filtered.bam
samtools sort POG_mapped_RNA_dud.filtered.bam > POG_mapped_RNA_dud.sorted.filtered.bam
samtools index POG_mapped_RNA_dud.sorted.filtered.bam

isoseq collapse --do-not-collapse-extra-5exons POG_mapped_RNA_dud.sorted.filtered.bam collapsed_isoseq_dud_reads.gff
