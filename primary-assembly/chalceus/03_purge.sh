#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name purge_dups
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH -o purgedups.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis
module load matplotlib/3.7.0-gfbf-2022b

# Working directory: /scratch/leuven/357/vsc35707/chalceus/
set -euo pipefail # will exit the script if something goes wrong

# ---------- INPUT ----------
PRI_ASM="Pogonus_T2T.asm.hic.p_ctg.fa"
ALT_ASM="Pogonus_T2T.asm.hic.a_ctg.fa"
PB_READS="pacbio/GC157810.fasta"
BIN_DIR="bin"  # path to purge_dups binaries
CPUs=36

# Check if bin directory exists
[ -d $BIN_DIR ] || mkdir -p $BIN_DIR

# ---------- STAGE 1: Purge from primary assembly ----------
echo "### STAGE 1: Purge_dups on primary assembly ###"

# Step 1a: Align PacBio reads to primary assembly
minimap2 -t $CPUs -xasm20 $PRI_ASM $PB_READS | gzip -c - > primary.paf.gz

# Step 1b: Compute coverage stats
$BIN_DIR/pbcstat primary.paf.gz

# Step 1c: Calculate cutoffs
$BIN_DIR/calcuts PB.stat > cutoffs_primary 2> calcuts_primary.log

# Step 1d: Split primary assembly
$BIN_DIR/split_fa $PRI_ASM > ${PRI_ASM}.split

# Step 1e: Self-alignment
minimap2 -t $CPUs -xasm5 -DP ${PRI_ASM}.split ${PRI_ASM}.split | gzip -c - > ${PRI_ASM}.split.self.paf.gz

# Step 2: Purge
$BIN_DIR/purge_dups -2 -T cutoffs_primary -c PB.base.cov ${PRI_ASM}.split.self.paf.gz > dups_primary.bed 2> purge_dups_primary.log

# Step 3: Extract sequences
$BIN_DIR/get_seqs -e dups_primary.bed $PRI_ASM

# This creates:
# - purged.fa (purged primary contigs)
# - hap.fa (candidate haplotigs from primary)

# Step 4: Merge hap.fa with alternative assembly
cat hap.fa $ALT_ASM > merged_hap.fa

# ---------- STAGE 2: Refine merged haplotigs ----------
echo "### STAGE 2: Purge_dups on merged haplotigs ###"

# Step 1a: Align PacBio reads
minimap2 -t $CPUs -xasm20 merged_hap.fa $PB_READS | gzip -c - > merged_hap.paf.gz

# Step 1b: Compute stats
$BIN_DIR/pbcstat merged_hap.paf.gz

# Step 1c: Calculate cutoffs
$BIN_DIR/calcuts PB.stat > cutoffs_merged 2> calcuts_merged.log

# Step 1d: Split merged hap assembly
$BIN_DIR/split_fa merged_hap.fa > merged_hap.fa.split

# Step 1e: Self-alignment
minimap2 -t $CPUs -xasm5 -DP merged_hap.fa.split merged_hap.fa.split | gzip -c - > merged_hap.fa.split.self.paf.gz

# Step 2: Purge again
$BIN_DIR/purge_dups -2 -T cutoffs_merged -c PB.base.cov merged_hap.fa.split.self.paf.gz > dups_merged.bed 2> purge_dups_merged.log

# Step 3: Extract final haplotigs
$BIN_DIR/get_seqs -e dups_merged.bed merged_hap.fa

# Rename output for clarity
mv purged.fa purged_merged.fa
mv hap.fa hap_merged.fa

echo "### DONE: Final haplotig assembly is in purged_merged.fa ###"
