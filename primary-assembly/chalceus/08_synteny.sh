#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name synteny
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00
#SBATCH -o synteny.%j.out
#SBATCH -A lp_svbelleghem

awk '/^>/ { if (count++ >= 20) exit } { print }' P_chalceus_REF1_sorted.fa  > 20_P_chalceus_REF1.fa 

# variables
REFERENCE=20_P_chalceus_REF1.fa
SCAFFOLDS=11_ragtag.scaffold.sorted.fa
PREFIX=P_chal_SW_RagTag

# index files
REFERENCE_FAI=${REFERENCE}.fai
SCAFFOLDS_FAI=${SCAFFOLDS}.fai

# Activate environment
conda activate ntsynt

# Step 1: Run ntSynt
ntSynt $REFERENCE $SCAFFOLDS -p $PREFIX -t 36 -d 30

# Step 2: Block statistics
python denovo_synteny_block_stats.py --tsv $PREFIX.synteny_blocks.tsv --fai $REFERENCE_FAI $SCAFFOLDS_FAI

# Step 3: Sort blocks
python sort_ntsynt_blocks.py \
  --synteny_blocks $PREFIX.synteny_blocks.tsv \
  --sort_order $REFERENCE_FAI $SCAFFOLDS_FAI --fais > $PREFIX.synteny_blocks.sorted.tsv

# Step 4: Format for gggenomes
python format_blocks_gggenomes.py \
  --fai $REFERENCE_FAI $SCAFFOLDS_FAI \
  --prefix $PREFIX \
  --blocks $PREFIX.synteny_blocks.sorted.tsv \
  --length 100 \
  --colour $REFERENCE

# Step 5: Copy results to VSC
cp $PREFIX.links.tsv $VSC_DATA
cp $PREFIX.sequence_lengths.tsv $VSC_DATA


# to run locally 
# Rscript plot_synteny_blocks_gggenomes.R -s P_chal_SW.sequence_lengths.tsv -l P_chal_SW.links.tsv --scale 25000000 --p P_chal_SW
