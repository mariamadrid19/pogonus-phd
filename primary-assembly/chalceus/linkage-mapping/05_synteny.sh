#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name synteny
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00
#SBATCH -o synteny.%j.out
#SBATCH -A lp_svbelleghem

# variables
REFERENCE=Pterostichus_genome.fna
GENOME1=final.fasta
GENOME2=26_Calosoma_wilkesii_sorted.fna
PREFIX=carabids

# index files
REFERENCE_FAI=${REFERENCE}.fai
GENOME1_FAI=${GENOME1}.fai
GENOME2_FAI=${GENOME2}.fai

# Activate environment
conda activate ntsynt

# Step 1: Run ntSynt
ntSynt $REFERENCE $GENOME1 $GENOME2 -p $PREFIX -t 36 -d 15

# Step 2: Block statistics
python denovo_synteny_block_stats.py --tsv $PREFIX.synteny_blocks.tsv --fai $REFERENCE_FAI $GENOME1_FAI $GENOME2_FAI

# Step 3: Sort blocks
python sort_ntsynt_blocks.py \
  --synteny_blocks $PREFIX.synteny_blocks.tsv \
  --sort_order $REFERENCE_FAI $GENOME1_FAI $GENOME2_FAI --fais > $PREFIX.synteny_blocks.sorted.tsv

# Step 4: Format for gggenomes
python format_blocks_gggenomes.py \
  --fai $REFERENCE_FAI $GENOME1_FAI $GENOME2_FAI \
  --prefix $PREFIX \
  --blocks $PREFIX.synteny_blocks.sorted.tsv \
  --length 100 \
  --colour $REFERENCE

# Step 5: Copy results to VSC
cp $PREFIX.links.tsv $VSC_DATA
cp $PREFIX.sequence_lengths.tsv $VSC_DATA
