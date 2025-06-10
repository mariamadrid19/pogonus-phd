#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=blobtoolkit
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=150G
#SBATCH --time=72:00:00
#SBATCH -o blobtoolkit.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/decontamination/

# === Load environment and modules ===
conda activate decon_env
module load SAMtools/1.16.1-GCC-11.3.0

# === Set variables ===
THREADS=36

# File locations
OUTPUT=./decontam_output_regions
ASSEMBLY=$OUTPUT/assembly.cleaned.fasta

PACBIO=/scratch/leuven/357/vsc35707/chalceus/pacbio/bc2041.fastq.gz
NANOPORE=/scratch/leuven/357/vsc35707/chalceus/nanopore/GC157812.fastq.gz
KRAKEN_OUT=$OUTPUT/kraken2.out
TAXDUMP=./taxdump            # Should contain nodes.dmp and names.dmp

# BlobToolKit directories
BLOB_TAX=$OUTPUT/blobtools_tax
BLOB_DATASET=$OUTPUT/blobtools_dataset

mkdir -p $BLOB_TAX $BLOB_DATASET

# === Step 1: Index the cleaned assembly for mapping ===
echo "Indexing assembly..."
minimap2 -d $OUTPUT/assembly.cleaned.mmi $ASSEMBLY

# === Step 2: Map reads to the cleaned assembly ===

# PacBio mapping
echo "Mapping PacBio reads..."
minimap2 -t $THREADS -ax map-pb $OUTPUT/assembly.cleaned.mmi $PACBIO | \
  samtools view -@ $THREADS -Sb - | \
  samtools sort -@ $THREADS -o $OUTPUT/pacbio.bam
samtools index $OUTPUT/pacbio.bam

# Nanopore mapping
echo "Mapping Nanopore reads..."
minimap2 -t $THREADS -ax map-ont $OUTPUT/assembly.cleaned.mmi $NANOPORE | \
  samtools view -@ $THREADS -Sb - | \
  samtools sort -@ $THREADS -o $OUTPUT/nanopore.bam
samtools index $OUTPUT/nanopore.bam

# Merge all BAMs for coverage
echo "Merging all alignments..."
samtools merge -@ $THREADS $OUTPUT/all_reads.bam $OUTPUT/pacbio.bam $OUTPUT/nanopore.bam
samtools index $OUTPUT/all_reads.bam

# === Step 3: Add taxonomy from Kraken2 ===
echo "Adding taxonomic labels..."
blobtools2 addtax \
  --taxid-files $KRAKEN_OUT \
  --taxdump $TAXDUMP \
  --output $BLOB_TAX

# === Step 4: Create BlobToolKit dataset ===
echo "Creating BlobToolKit dataset..."
blobtools2 create \
  --fasta $ASSEMBLY \
  --bam $OUTPUT/all_reads.bam \
  --taxid-files $BLOB_TAX/taxid.json \
  --taxdump $TAXDUMP \
  --output $BLOB_DATASET

# === Step 5: Generate static HTML plots ===
echo "Generating BlobToolKit plots..."
blobtools2 plot --directory $BLOB_DATASET

echo "=== BlobToolKit analysis complete ==="
