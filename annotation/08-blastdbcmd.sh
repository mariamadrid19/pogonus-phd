#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name blastdbcmd
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128G
#SBATCH --time=72:00:00
#SBATCH --output=logs/blastdbcmd_%A_%a.out
#SBATCH --error=logs/blastdbcmd_%A_%a.err
#SBATCH --array=0-17
#SBATCH -A lp_svbelleghem

module load BLAST+/2.13.0-gompi-2022a

export BLASTDB=/scratch/leuven/357/vsc35707/eviann/more-data/nr

mkdir -p blast_out

# Get list of input files into an array
FILES=(insecta_part_*)
INPUT_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

# Output name based on input file
OUT_FILE="blast_out/${INPUT_FILE}.fasta"

# Run blastdbcmd
blastdbcmd -db ./nr -entry_batch "$INPUT_FILE" -out "$OUT_FILE" -outfmt %f
