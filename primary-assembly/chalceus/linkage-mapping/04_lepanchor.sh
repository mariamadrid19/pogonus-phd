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

# Filter and process map.clean
awk '($5>=1.0)' map.clean | awk -f $LEPANCHOR/cleanmap.awk | awk -f $LEPANCHOR/joinIntervals.awk $GENOME.sizes - > map.bed

# Create .bed files for each chromosome 
awk '{fn="chr" $4 ".bed"; print $1"\t"$2"\t"$3 > fn}' map.bed

# Run OrderMarkers2 on chromosomes 1 to 11
for X in {1..11}; do
  java -cp $LEPMAP/bin OrderMarkers2 map=map5.txt data=data_f.call chromosome=$X recombination1=0 > order${X}.txt
done

# Prepare input for PlaceAndOrientContigs
for x in {1..11}; do
  awk -vn=$x '(NR==FNR){map[NR-1]=$0} (NR!=FNR && /^[^#]/){print map[$1],n,$2,$3}' snps.txt order$x.txt > order$x.m.input
done

# Run PlaceAndOrientContigs
for X in {1..11}; do
  java -cp $LEPANCHOR/bin/ PlaceAndOrientContigs map=order$X.m.input bed=chr$X.bed noIntervals=1 > chr$X.la 2> chr$X.la.err
done

# Generate .agp files
for f in chr*.la; do
  chrnum=$(echo $f | sed 's/[^0-9]//g')
  awk -vlg=$chrnum -f $LEPANCHOR/makeagp_full2.awk $f > chr$chrnum.agp
done

# Create the final AGP file
cat chr*.agp > final.agp

# Create the final fasta file
awk -f $LEPANCHOR/makefasta.awk $GENOME chr*.agp > final.fasta

module load SAMtools/1.16.1-GCC-11.3.0

# Index the fasta file
samtools faidx final.fasta
