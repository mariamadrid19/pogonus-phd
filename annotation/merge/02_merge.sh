#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=merge_LW_annotations
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=72G
#SBATCH -A lp_svbelleghem
#SBATCH -o merge_LW_annotations.%j.out
#SBATCH -e merge_LW_annotations.%j.err

set -euo pipefail

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate agat-merge

output_dir=/scratch/leuven/357/vsc35707/annotation/LW-annotation/merge-annotations

braker="${output_dir}/braker.normalized.gff3"
eviann="${output_dir}/eviann.normalized.gff3"

merged="${output_dir}/Pchalceus_LW.BRAKER_EviAnn.merged.gff3"
merged_stats="${output_dir}/Pchalceus_LW.BRAKER_EviAnn.merged.statistics.txt"

mkdir -p "$output_dir"
cd "$output_dir"

# Verify software.
perl -MFile::ShareDir -e \
    'print "File::ShareDir loaded successfully\n"'

command -v agat_sp_complement_annotations.pl
command -v agat_sp_statistics.pl

# Verify inputs.
test -s "$braker"
test -s "$eviann"

echo "BRAKER input: $braker"
echo "EviAnn input: $eviann"

# Count genes before merging.
braker_genes=$(
    awk -F '\t' '
        $0 !~ /^#/ && tolower($3)=="gene" {n++}
        END {print n+0}
    ' "$braker"
)

eviann_genes=$(
    awk -F '\t' '
        $0 !~ /^#/ && tolower($3)=="gene" {n++}
        END {print n+0}
    ' "$eviann"
)

echo "BRAKER genes: $braker_genes"
echo "EviAnn genes: $eviann_genes"

# Preserve BRAKER wherever annotations conflict.
# Add EviAnn genes only when allowed by AGAT's complement rules.
agat_sp_complement_annotations.pl \
    --ref "$braker" \
    --add "$eviann" \
    --out "$merged"

test -s "$merged"

# Generate statistics for the merged annotation.
agat_sp_statistics.pl \
    --gff "$merged" \
    --output "$merged_stats"

merged_genes=$(
    awk -F '\t' '
        $0 !~ /^#/ && tolower($3)=="gene" {n++}
        END {print n+0}
    ' "$merged"
)

added_genes=$((merged_genes - braker_genes))

echo "Merged genes: $merged_genes"
echo "Genes added relative to BRAKER: $added_genes"
echo "Merged annotation: $merged"
echo "Statistics: $merged_stats"
