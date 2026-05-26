#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=lepanchor_clean
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1:00:00
#SBATCH -o lepanchor_clean.%j.out
#SBATCH -A lp_edu_eeg_2026

module load Java/25.36
conda activate variant_tools

set -euo pipefail

cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap/results-nonrepeat

LEPANCHOR="/data/leuven/357/vsc35707/LepAnchor"

MAP_PREFIX="map12.nonrepeat"
NCHR=11
MINPUT_SUFFIX="_m"

# New thresholds
MIN_MARKERS=30
MIN_SPAN_BP=10000000
MIN_SPAN_FRAC=0.50

MAP_BED="${MAP_PREFIX}.bed"

# overwrite previous cleaned outputs
rm -rf clean_vertical_reports
rm -rf clean_las
rm -rf clean_agp
rm -f ${MAP_PREFIX}_chr*_m.clean.input

mkdir -p logs clean_vertical_reports clean_las clean_agp

# =========================================================
# Remove suspicious vertical-line marker groups
# =========================================================

for X in $(seq 1 "$NCHR"); do
  IN="${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.input"
  OUT="${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.clean.input"
  BAD="clean_vertical_reports/${MAP_PREFIX}_chr${X}.bad_vertical_markers.txt"
  REPORT="clean_vertical_reports/${MAP_PREFIX}_chr${X}.vertical_groups_report.txt"

  echo "Cleaning chromosome ${X}: ${IN}"

  awk -v min_n="$MIN_MARKERS" \
      -v min_span_bp="$MIN_SPAN_BP" \
      -v min_span_frac="$MIN_SPAN_FRAC" \
      -v report="$REPORT" \
      -v badfile="$BAD" \
      -v outfile="$OUT" '
    {
      scaffold=$1
      bp=$2 + 0
      lg=$3
      female=$5

      key=scaffold OFS lg OFS female

      n[key]++

      if (!(key in minbp) || bp < minbp[key]) minbp[key]=bp
      if (!(key in maxbp) || bp > maxbp[key]) maxbp[key]=bp

      if (!(scaffold in scaffold_min) || bp < scaffold_min[scaffold]) scaffold_min[scaffold]=bp
      if (!(scaffold in scaffold_max) || bp > scaffold_max[scaffold]) scaffold_max[scaffold]=bp

      line[NR]=$0
      linekey[NR]=key
      lineid[NR]=scaffold OFS bp
    }

    END {
      print "scaffold", "LG", "female_cM", "n_markers", "min_bp", "max_bp", "span_bp", "scaffold_marker_span", "span_fraction" > report

      for (k in n) {
        split(k, a, OFS)
        scaffold=a[1]

        span=maxbp[k]-minbp[k]
        scaffold_span=scaffold_max[scaffold]-scaffold_min[scaffold]
        frac=(scaffold_span > 0 ? span/scaffold_span : 0)

        if (n[k] >= min_n && span >= min_span_bp && frac >= min_span_frac) {
          bad[k]=1
          print k, n[k], minbp[k], maxbp[k], span, scaffold_span, frac >> report
        }
      }

      for (i=1; i<=NR; i++) {
        if (linekey[i] in bad) {
          print lineid[i] > badfile
        } else {
          print line[i] > outfile
        }
      }
    }
  ' "$IN"

  echo "  removed markers: $(wc -l < "$BAD")"
done

# =========================================================
# Run PlaceAndOrientContigs
# =========================================================

for X in $(seq 1 "$NCHR"); do
  INPUT="${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.clean.input"

  echo "Running PlaceAndOrientContigs for chromosome ${X}"

  java -cp "$LEPANCHOR/bin/" PlaceAndOrientContigs \
    map="$INPUT" \
    bed="$MAP_BED" \
    chromosome="$X" \
    noIntervals=1 \
    > "clean_las/${MAP_PREFIX}_chr${X}.clean.la" \
    2> "clean_las/${MAP_PREFIX}_chr${X}.clean.la.err"
done

# =========================================================
# Generate cleaned AGP files
# =========================================================

for X in $(seq 1 "$NCHR"); do
  echo "Generating AGP for chromosome ${X}"

  awk -vlg="$X" -f "$LEPANCHOR/makeagp_full2.awk" \
    "clean_las/${MAP_PREFIX}_chr${X}.clean.la" \
    > "clean_agp/${MAP_PREFIX}_chr${X}.clean.agp"
done

echo "Done."
echo "Thresholds used:"
echo "  MIN_MARKERS=${MIN_MARKERS}"
echo "  MIN_SPAN_BP=${MIN_SPAN_BP}"
echo "  MIN_SPAN_FRAC=${MIN_SPAN_FRAC}"
echo ""
echo "Cleaned input files: ${MAP_PREFIX}_chr*_m.clean.input"
echo "Reports: clean_vertical_reports/"
echo "LepAnchor outputs: clean_las/"
echo "AGP files: clean_agp/"
