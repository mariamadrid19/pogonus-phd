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

# GC157812 (barcode 31) = 58 Gb (short-winged, SW)
# GC157813 (barcode 39) = 20 Gb (long-winged, LW -> probably P. gilvipes instead of P. chalceus like we expected) 

# Count reads in original files
zcat fastq_pass/barcode31/*.fastq.gz | echo $((`wc -l` / 4))
zcat fastq_pass/barcode39/*.fastq.gz | echo $((`wc -l` / 4))

# Count reads in merged file
zcat GC157812.fastq.gz | echo $((`wc -l` / 4))
zcat GC157813.fastq.gz | echo $((`wc -l` / 4))

zcat GC157812.fastq.gz | awk '{c++} END{print "Total lines:", c, "Reads:", c/4}'
# total of 13533617 reads, 54134468 lines
zcat GC157813.fastq.gz | awk '{c++} END{print "Total lines:", c, "Reads:", c/4}'
# 5913014 reads, 23652056 lines

module load FastQC/0.12.1-Java-11
# check quality of reads with fastqc, not the most realiable since fastqc is made for short reads
fastqc GC157812.fastq.gz -t 32
fastqc GC157813.fastq.gz -t 32

# try with a qc designed for ONT reads
module load NanoPlot/1.42.0-foss-2022a
NanoPlot --fastq GC157812.fastq.gz -o nanoplot_SW
NanoPlot --fastq GC157813.fastq.gz -o nanoplot_GILV

# can also try with pycoQC
conda activate pycoQC
pycoQC -f sequencing_summary_PAW02052_066d7b83_dd7ec26e.txt -o pycoqc_report_dd7ec26e.html
pycoQC -f sequencing_summary_PAW02052_a0329603_c7340ae4.txt -o pycoqc_report_c7340ae4.html
