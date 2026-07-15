#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=pb_qc
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH --mem=120G
#SBATCH -o pb_qc.%j.out
#SBATCH -e pb_qc.%j.err
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/sw-assembly

source "$(conda info --base)/etc/profile.d/conda.sh"

mkdir -p fastq fastqc nanoplot jellyfish

# Convert both barcode BAM files separately
bam2fastq \
    -o fastq/bc2041 \
    hifi_reads/m84247_241127_161129_s3.hifi_reads.bc2041.bam

bam2fastq \
    -o fastq/bc2042 \
    hifi_reads/m84247_241127_161129_s3.hifi_reads.bc2042.bam

# Combine both barcode FASTQs because they come from the same individual
cat \
    fastq/bc2041.fastq.gz \
    fastq/bc2042.fastq.gz \
    > fastq/P_chalceus_HiFi_merged.fastq.gz

# Check that the merged FASTQ is readable
gzip -t fastq/P_chalceus_HiFi_merged.fastq.gz

module load FastQC/0.12.1-Java-11

fastqc \
    fastq/P_chalceus_HiFi_merged.fastq.gz \
    -t 32 \
    -o fastqc

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate thesis

module load NanoPlot/1.42.0-foss-2022a

NanoPlot \
    --fastq fastq/P_chalceus_HiFi_merged.fastq.gz \
    --threads 32 \
    -o nanoplot/P_chalceus_HiFi_merged

# Convert merged FASTQ to FASTA
seqtk seq -a \
    fastq/P_chalceus_HiFi_merged.fastq.gz \
    > P_chalceus_HiFi_merged.fasta

module purge
conda activate thesis

# Count canonical 31-mers
jellyfish count \
    -m 31 \
    -s 10G \
    -t 32 \
    -C \
    -o jellyfish/P_chalceus_k31.jf \
    P_chalceus_HiFi_merged.fasta

# Generate histogram for GenomeScope2
jellyfish histo \
    -t 32 \
    jellyfish/P_chalceus_k31.jf \
    > jellyfish/P_chalceus_k31.histo

echo "Pipeline completed successfully."
