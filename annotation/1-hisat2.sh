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

indID=$((SLURM_ARRAY_TASK_ID -1))

ind=( )

hisat2 --rna-strandness RF -x sorted_prim_dud -1 /RawData/$(echo "${ind[indID]}")_1.fq.gz -2 /RawData/$(echo "${ind[indID]}")_2.fq.gz -S $(echo "${ind[indID]}").sam
 
samtools view -h -q 30 -b $(echo "${ind[indID]}").sam > $(echo "${ind[indID]}")._filtered.bam
samtools sort $(echo "${ind[indID]}")._filtered.bam > $(echo "${ind[indID]}").filtered.bam
samtools index $(echo "${ind[indID]}").filtered.bam
 
