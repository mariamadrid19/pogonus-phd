#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name purge_dups
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH --output=logs/purgedups.%j.out
#SBATCH --error=logs/purgedups.%j.err
#SBATCH -A lp_svbelleghem

conda activate thesis
module load matplotlib/3.7.0-gfbf-2022b

# Working directory: /scratch/leuven/357/vsc35707/olivaceus/
set -euo pipefail 
# will exit the script if something goes wrong

# ---------- INPUT ----------
PRI_ASM="Pogonus_olivaceus.asm.bp.p_ctg.fa"
PB_READS="P_olivaceus_HiFi_top30.fasta"
CPUs=36

# ---------- STAGE 1: Purge from primary assembly ----------
echo "### STAGE 1: Purge_dups on primary assembly ###"

# Step 1a: Align PacBio reads to primary assembly
minimap2 -t $CPUs -xasm20 $PRI_ASM $PB_READS | gzip -c - > primary.paf.gz

# Step 1b: Compute coverage stats
pbcstat primary.paf.gz

# Step 1c: Calculate cutoffs
calcuts PB.stat > cutoffs_primary 2> calcuts_primary.log

# Step 1d: Split primary assembly
split_fa $PRI_ASM > ${PRI_ASM}.split

# Step 1e: Self-alignment
minimap2 -t $CPUs -xasm5 -DP ${PRI_ASM}.split ${PRI_ASM}.split | gzip -c - > ${PRI_ASM}.split.self.paf.gz

# Step 2: Purge
purge_dups -2 -T cutoffs_primary -c PB.base.cov ${PRI_ASM}.split.self.paf.gz > dups_primary.bed 2> purge_dups_primary.log

# Step 3: Extract sequences
get_seqs -e dups_primary.bed $PRI_ASM

# ---------- STAGE 2: Check completeness ----------
