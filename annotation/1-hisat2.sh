#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name hisat2
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH --output=hisat2.%j.out
#SBATCH -A lp_svbelleghem

module load HISAT2/2.2.1-gompi-2021b
module load SAMtools

hisat2-build sorted_prim_dud.fasta sorted_prim_dud

hisat2 -f -x sorted_prim_dud -U POG_IsoSeq_HiFi_demux.fastq -S Dud_mapped_IsoSeq.sam

samtools view -h -q 30 -b Dud_mapped_IsoSeq.sam > Dud_mapped_IsoSeq_filtered.bam
samtools sort Dud_mapped_IsoSeq_filtered.bam > Dud_mapped_IsoSeq.filtered.bam
samtools index Dud_mapped_IsoSeq.filtered.bam
