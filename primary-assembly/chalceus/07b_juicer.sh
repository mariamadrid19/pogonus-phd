#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=juicer
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00
#SBATCH -o juicer.%j.out
#SBATCH -A lp_svbelleghem

MAPPED='mapped2.PT.bam'
AGP='yahs.out_scaffolds_final.agp'
INDEX_CONTIGS='purged.fa.fai'

# Uses the mapped HiC to the scaffolds, the YaHS AGP output file, and the index file of the purged contigs to generate the file needed for juicer
juicer pre -a -o out_JBAT pre $MAPPED $AGP $INDEX_CONTIGS > out_JBAT.log 2 > &1

# Generates a HiC contact map file that can be loaded by Juicebox for manual editing
(java -Xmx48000m  -Djava.awt.headless=true -jar juicer_tools.2.20.00.jar pre out_JBAT.txt out_JBAT.hic.part <(cat out_JBAT.log  | grep PRE_C_SIZE | awk '{print $2" "$3}')) && (mv out_JBAT.hic.part out_JBAT.hic)
