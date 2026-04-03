#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=lepanchor
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH -o logs/lepanchor.%j.out
#SBATCH -A lp_svbelleghem

module load Java/21.0.8
conda activate variant_tools

set -euo pipefail

# Change to the working directory
cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap

# Define directories
LEPANCHOR="/data/leuven/357/vsc35707/LepAnchor"
LEPMAP="/data/leuven/357/vsc35707/LepMap3"
GENOME="/scratch/leuven/357/vsc35707/linkage-mapping/genome/Pchalceus_LW.fa"

# =========================
# User-defined variables
# =========================
MAP_PREFIX="map12"     # base name, without .txt
NCHR=11
USE_MRECOM0=1              # 1 = use *_mrecom0.txt files, 0 = use default OrderMarkers2 files

MAP_FILE="${MAP_PREFIX}.txt"
CLEANMAP_INPUT="${MAP_PREFIX}.cleanMap.input"
CLEANMAP_SORTED="${MAP_PREFIX}.cleanMap.sorted.input"
MAP_CLEAN="${MAP_PREFIX}.clean"
MAP_BED="${MAP_PREFIX}.bed"

# Suffix for OrderMarkers2 files
if [[ "$USE_MRECOM0" -eq 1 ]]; then
    ORDER_SUFFIX="_mrecom0"
    MINPUT_SUFFIX="_m"
else
    ORDER_SUFFIX=""
    MINPUT_SUFFIX=""
fi

echo "Using MAP_PREFIX=${MAP_PREFIX}"
echo "Using MAP_FILE=${MAP_FILE}"
echo "Using USE_MRECOM0=${USE_MRECOM0}"

# How many markers per LG?
cut -f1 "$MAP_FILE" | sort -n | uniq -c

# Create map input for CleanMap
paste snps.txt "$MAP_FILE" | awk '(NR>1)' > "$CLEANMAP_INPUT"

# cleanMap.input should be sorted by contig and position for CleanMap
sort -V -k1,1 -k2,2n "$CLEANMAP_INPUT" > "$CLEANMAP_SORTED"

# Clean map file
java -cp "$LEPANCHOR/bin/" CleanMap map="$CLEANMAP_SORTED" > "$MAP_CLEAN"

# Generate genome sizes file if needed
awk -f "$LEPANCHOR/contigLength.awk" "$GENOME" > "$GENOME.sizes"

# Generate .bed file for the entire genome
java -cp "$LEPANCHOR/bin/" Map2Bed map="$MAP_CLEAN" contigLength="$GENOME.sizes" > "$MAP_BED"

# Prepare input for PlaceAndOrientContigs from OrderMarkers2 output
for X in $(seq 1 "$NCHR"); do
  awk -vn="$X" '
    (NR==FNR){map[NR-1]=$0}
    (NR!=FNR && /^[^#]/){print map[$1], n, $2, $3}
  ' snps.txt "${MAP_PREFIX}_chr${X}${ORDER_SUFFIX}.txt" > "${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.input"
done

# =========================
# Run PlaceAndOrientContigs
# =========================
for X in $(seq 1 "$NCHR"); do
  java -cp "$LEPANCHOR/bin/" PlaceAndOrientContigs \
    map="${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.input" \
    bed="$MAP_BED" \
    chromosome="$X" \
    noIntervals=1 \
    > "${MAP_PREFIX}_chr${X}.la" 2> "${MAP_PREFIX}_chr${X}.la.err"
done

# =========================
# Generate .agp files
# =========================
for X in $(seq 1 "$NCHR"); do
  awk -vlg="$X" -f "$LEPANCHOR/makeagp_full2.awk" \
    "${MAP_PREFIX}_chr${X}.la" > "${MAP_PREFIX}_chr${X}.agp"
done
