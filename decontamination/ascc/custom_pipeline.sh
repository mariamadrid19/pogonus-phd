#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --partition=batch_long
#SBATCH --job-name=decontamination
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=150G
#SBATCH --time=100:00:00
#SBATCH -o decon.%j.out
#SBATCH -A lp_svbelleghem

# === Activate conda + modules ===
conda activate decon_env
module load BLAST+/2.13.0-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0

cd /scratch/leuven/357/vsc35707/decontamination/

# === Variables ===
THREADS=36
ASSEMBLY=purged.fa
KRAKEN_DB=./kraken2db
NR_DB=./nr_dmnd
TAXDUMP=./taxdump
OUTPUT=decontam_output_regions
TARGET_TAXID=235516  # Pogonus chalceus

mkdir -p $OUTPUT

# === Step 1: Kraken2 classification ===
kraken2 --db $KRAKEN_DB --threads $THREADS \
  --report $OUTPUT/kraken2.report \
  --output $OUTPUT/kraken2.out \
  $ASSEMBLY

# Extract contigs assigned to non-target taxIDs (contaminants)
awk -v target=$TARGET_TAXID '$6 != target {print $2}' $OUTPUT/kraken2.out > $OUTPUT/kraken_contaminants.txt

# === Step 2: Tiara classification ===
tiara -i $ASSEMBLY -o $OUTPUT/tiara.tsv -t $THREADS
awk '$2 != "eukaryote"' $OUTPUT/tiara.tsv | cut -f1 > $OUTPUT/tiara_contaminants.txt

# === Step 3: DIAMOND BLASTX (nr) ===
diamond blastx -d $NR_DB -q $ASSEMBLY -o $OUTPUT/diamond.tsv \
  -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  -k 1 --evalue 1e-5 --threads $THREADS
