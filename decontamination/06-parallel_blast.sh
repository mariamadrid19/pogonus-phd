#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name blast_scafs
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=32
#SBATCH --time=04:00:00
#SBATCH --output=blast_scafs.%j.out
#SBATCH --array=1-150
#SBATCH -A lp_svbelleghem

# Input files
ASSEMBLY="sorted_prim_dud.fasta"
SCAFFOLDS="scaffold_names.txt"

# Extract the specific scaffold name based on SLURM array ID
SCAF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SCAFFOLDS")

echo "Processing scaffold: $SCAF"
echo "Assembly file: $ASSEMBLY"

# Extract the scaffold sequence into a FASTA file
awk -v scaffold=">${SCAF}" '
BEGIN {print_seq=0}
$0 ~ scaffold {print_seq=1; print $0; next}
$0 ~ /^>/ {print_seq=0}
print_seq' "$ASSEMBLY" > "${SCAF}.fasta"

module load BLAST+/2.13.0-gompi-2022a

export BLASTDB=/scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud/nt

# Run BLAST for this scaffold
blastn -db nt \
       -query "${SCAF}.fasta" \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 32 \
       -taxidlist microbial.ids \
       -out "${SCAF}.ncbi.blastn.out"

# Clean up scaffold FASTA file
rm "${SCAF}.fasta"

# Directory containing BLAST output files
OUTPUT_DIR="/scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud"
FINAL_OUTPUT="sorted_prim_dud.ncbi.blastn.out"

# Combine all *.ncbi.blastn.out files
cat ${OUTPUT_DIR}/*.ncbi.blastn.out > ${OUTPUT_DIR}/${FINAL_OUTPUT}

echo "All BLAST outputs combined into ${FINAL_OUTPUT}"
