#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=ont_cat
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=32
#SBATCH --mem-per-cpu=8G 
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o ont_cat.%j.out

# concatenate and compress all the fastq files that "pass" quality checks, ONT reads
zcat fastq_pass/barcode31/*.fastq.gz | pigz -p 12 > GC157812.fastq.gz
zcat fastq_pass/barcode39/*.fastq.gz | pigz -p 12 > GC157813.fastq.gz

# Count reads in original files
zcat fastq_pass/barcode31/*.fastq.gz | echo $((`wc -l` / 4))
zcat fastq_pass/barcode39/*.fastq.gz | echo $((`wc -l` / 4))

# Count reads in merged file
zcat GC157812.fastq.gz | echo $((`wc -l` / 4))
zcat GC157813.fastq.gz | echo $((`wc -l` / 4))

zcat GC157812.fastq.gz | awk '{c++} END{print "Total lines:", c, "Reads:", c/4}'
zcat GC157813.fastq.gz | awk '{c++} END{print "Total lines:", c, "Reads:", c/4}'

module load FastQC/0.12.1-Java-11
# check quality of reads with fastqc
fastqc GC157812.fastq.gz -t 32
fastqc GC157813.fastq.gz -t 32
