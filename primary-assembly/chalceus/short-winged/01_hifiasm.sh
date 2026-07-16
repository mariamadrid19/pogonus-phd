#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --partition=bigmem
#SBATCH --job-name=Pog_hifiasm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --cpus-per-task=1
#SBATCH --mem=1200G
#SBATCH --time=72:00:00
#SBATCH --output=logs/Pog_hifiasm.%j.out
#SBATCH --error=logs/Pog_hifiasm.%j.err
#SBATCH --account=lp_edu_eeg_2026

# Working directory
WORKDIR="/scratch/leuven/357/vsc35707/sw-assembly"
cd "$WORKDIR"

# Input files
HIFI="${WORKDIR}/P_chalceus_HiFi_merged.fasta"

HIC_R1="${WORKDIR}/hiC_reads/GC164167_TTCCAAGG-CCTTGTAG_S10_L001_R1_001.fastq.gz"
HIC_R2="${WORKDIR}/hiC_reads/GC164167_TTCCAAGG-CCTTGTAG_S10_L001_R2_001.fastq.gz"

# Output directory and prefix
OUTDIR="${WORKDIR}/hifiasm_HiFi_HiC"
PREFIX="${OUTDIR}/Pogonus_SW.asm"

mkdir -p "$OUTDIR"

# Activate the environment containing hifiasm
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate thesis

# Check that all required inputs exist and are non-empty
for file in "$HIFI" "$HIC_R1" "$HIC_R2"
do
    if [[ ! -s "$file" ]]; then
        echo "ERROR: Input file missing or empty: $file" >&2
        exit 1
    fi
done

echo "Hifiasm version:"
hifiasm --version

echo "Starting HiFi + Hi-C assembly"
date

hifiasm \
    -o "$PREFIX" \
    -t 24 \
    --h1 "$HIC_R1" \
    --h2 "$HIC_R2" \
    --n-hap 2 \
    "$HIFI" \
    2> "${PREFIX}.log"

echo "Hifiasm finished"
date

# Confirm that the expected primary assembly graph exists
PRIMARY_GFA="${PREFIX}.hic.p_ctg.gfa"
HAP1_GFA="${PREFIX}.hic.hap1.p_ctg.gfa"
HAP2_GFA="${PREFIX}.hic.hap2.p_ctg.gfa"

if [[ ! -s "$PRIMARY_GFA" ]]; then
    echo "ERROR: Expected primary GFA was not created: $PRIMARY_GFA" >&2
    exit 1
fi

# Convert primary assembly GFA to FASTA
awk '/^S/{print ">"$2; print $3}' \
    "$PRIMARY_GFA" \
    > "${PREFIX}.hic.p_ctg.fa"

# Convert haplotype assemblies if they were generated
if [[ -s "$HAP1_GFA" ]]; then
    awk '/^S/{print ">"$2; print $3}' \
        "$HAP1_GFA" \
        > "${PREFIX}.hic.hap1.p_ctg.fa"
fi

if [[ -s "$HAP2_GFA" ]]; then
    awk '/^S/{print ">"$2; print $3}' \
        "$HAP2_GFA" \
        > "${PREFIX}.hic.hap2.p_ctg.fa"
fi

# Basic assembly summaries
for fasta in \
    "${PREFIX}.hic.p_ctg.fa" \
    "${PREFIX}.hic.hap1.p_ctg.fa" \
    "${PREFIX}.hic.hap2.p_ctg.fa"
do
    if [[ -s "$fasta" ]]; then
        echo
        echo "Assembly: $fasta"

        awk '
            /^>/ {
                n++
                next
            }
            {
                total += length($0)
            }
            END {
                print "Contigs:", n
                print "Total bp:", total
            }
        ' "$fasta"
    fi
done

echo
echo "Assembly completed successfully."
echo "Primary assembly: ${PREFIX}.hic.p_ctg.fa"
echo "Haplotype 1:     ${PREFIX}.hic.hap1.p_ctg.fa"
echo "Haplotype 2:     ${PREFIX}.hic.hap2.p_ctg.fa"
