#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=inspect_eviann
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH -A lp_svbelleghem
#SBATCH -o inspect_eviann.%j.out
#SBATCH -e inspect_eviann.%j.err

set -euo pipefail

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate agat-merge

input_dir=/scratch/leuven/357/vsc35707/annotation/LW-annotation/
output_dir=/scratch/leuven/357/vsc35707/annotation/LW-annotation/merge-annotations

eviann_annotation="${input_dir}/eviann-LW/Pchalceus_LW_final.fasta.pseudo_label.gff"
braker_annotation="${input_dir}/functional/results/07_annotated_annotation/braker.annotated.gff3"

normalized_eviann_annotation="${output_dir}/eviann.normalized.gff3"
normalized_braker_annotation="${output_dir}/braker.normalized.gff3"

mkdir -p "$output_dir"
cd "$output_dir"

# Verify programs and Perl dependencies before starting.
perl -MFile::ShareDir -e \
    'print "File::ShareDir loaded successfully\n"'

command -v agat_convert_sp_gxf2gxf.pl
command -v agat_sp_statistics.pl
command -v agat_sp_compare_two_annotations.pl

# Verify inputs.
test -s "$braker_annotation"
test -s "$eviann_annotation"

# Step 1: Normalize both annotations.
agat_convert_sp_gxf2gxf.pl \
    --gff "$braker_annotation" \
    --output "$normalized_braker_annotation"

agat_convert_sp_gxf2gxf.pl \
    --gff "$eviann_annotation" \
    --output "$normalized_eviann_annotation"

# Step 2: Inspect feature counts.
agat_sp_statistics.pl \
    --gff "$normalized_braker_annotation" \
    --output "${output_dir}/braker.statistics.txt"

agat_sp_statistics.pl \
    --gff "$normalized_eviann_annotation" \
    --output "${output_dir}/eviann.statistics.txt"

# Step 3: Compare loci, including possible splits and fusions.
# Current AGAT uses --gff1 and --gff2.
agat_sp_compare_two_annotations.pl \
    --gff1 "$normalized_braker_annotation" \
    --gff2 "$normalized_eviann_annotation" \
    --output "${output_dir}/comparison"

# Step 4: Convert to GTF for transcript-level comparison.
agat_convert_sp_gff2gtf.pl \
    --gff "$normalized_braker_annotation" \
    --output "${output_dir}/braker.annotation.gtf"

agat_convert_sp_gff2gtf.pl \
    --gff "$normalized_eviann_annotation" \
    --output "${output_dir}/eviann.annotation.gtf"

# Step 5: Run gffcompare.
# Use its absolute path if it is not installed in the active environment.
gffcompare_bin=/data/leuven/357/vsc35707/gffcompare/gffcompare

test -x "$gffcompare_bin"

"$gffcompare_bin" \
    -r "${output_dir}/braker.annotation.gtf" \
    -o "${output_dir}/eviann-vs-braker-LW" \
    "${output_dir}/eviann.annotation.gtf"

echo "Comparison completed successfully"
