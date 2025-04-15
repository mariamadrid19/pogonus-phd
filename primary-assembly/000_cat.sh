#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=ont_cat
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=8G 
#SBATCH --time=6:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o ont_cat.%j.out

module load pigz/2.8-GCCcore-13.3.0

# concatenate and compress all the fastq files that "pass" quality checks, ONT reads
zcat fastq_pass/barcode31/*.fastq.gz | pigz -p $SLURM_CPUS_PER_TASK > GC157812.fastq.gz
zcat fastq_pass/barcode39/*.fastq.gz | pigz -p $SLURM_CPUS_PER_TASK > GC157813.fastq.gz
