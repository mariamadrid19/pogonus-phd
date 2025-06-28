#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name lepanchor
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o lepanchor.%j.out

# Load Java (needed to run LEPMAP and LEPANCHOR)
module load Java/21.0.2

# Change to the working directory
cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap

# Define directories
LEPANCHOR="/data/leuven/357/vsc35707/LepAnchor"
LEPMAP="/data/leuven/357/vsc35707/LepMap3"
GENOME="/scratch/leuven/357/vsc35707/linkage-mapping/genome/P_chalceus_broken.fa"

# How many markers per LG?
cut -f 1 map5.txt|sort -n|uniq -c

# Clean map file
java -cp $LEPANCHOR/bin/ CleanMap map=cleanMap.input > map.clean

# Generate .bed file for the entire genome
java -cp $LEPANCHOR/bin Map2Bed map=map.clean contigLength=$GENOME.sizes > map.bed

# Create .bed files for each chromosome 
for i in {1..12}; do
  awk -v chr=$i '$5 == chr {print $1 "\t" $2 "\t" $3}' map.bed > chr${i}.bed
done

# Run OrderMarkers2 on chromosomes 1 to 11
for X in {1..12}; do
  java -cp $LEPMAP/bin OrderMarkers2 map=map5.txt data=data_f.call chromosome=$X recombination1=0 > order${X}.txt
done

# Prepare input for PlaceAndOrientContigs
for X in {1..12}; do
  awk -vn=$X '(NR==FNR){map[NR-1]=$0} (NR!=FNR && /^[^#]/){print map[$1],n,$2,$3}' snps.txt order$X.txt > order$X.m.input
done

# Run PlaceAndOrientContigs
for X in {1..12}; do
  java -cp $LEPANCHOR/bin/ PlaceAndOrientContigs map=order$X.m.input bed=chr$X.bed noIntervals=1 > M_chr$X.la 2> M_chr$X.la.err
done

# Generate .agp files
for i in {1..12}; do
  awk -vlg=$i -f $LEPANCHOR/makeagp_full2 chr$i.la >chr$i.agp
done

# Create the final fasta file
awk -f $LEPANCHOR/makefasta.awk $GENOME chr*.agp > final.fasta

module load SAMtools/1.16.1-GCC-11.3.0

# Index the fasta file
samtools faidx final.fasta
