# =========================
# Remove suspicious vertical-line marker groups from all chromosomes
# Criteria:
# same scaffold + same LG + same female cM
# AND enough markers
# AND large physical span
# =========================

MIN_MARKERS=5
MIN_SPAN_BP=1000000
MIN_SPAN_FRAC=0.10

FILTERED_INPUT_LIST="filtered_minput_files.txt"
> "$FILTERED_INPUT_LIST"

for X in $(seq 1 "$NCHR"); do
  IN="${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.input"
  OUT="${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.clean.input"
  BAD="${MAP_PREFIX}_chr${X}.bad_vertical_markers.txt"
  REPORT="${MAP_PREFIX}_chr${X}.vertical_groups_report.txt"

  awk -v min_n="$MIN_MARKERS" \
      -v min_span_bp="$MIN_SPAN_BP" \
      -v min_span_frac="$MIN_SPAN_FRAC" '
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
  ' report="$REPORT" badfile="$BAD" outfile="$OUT" "$IN"

  echo "$OUT" >> "$FILTERED_INPUT_LIST"

  echo "Chromosome $X"
  echo "  input: $IN"
  echo "  output: $OUT"
  echo "  removed markers:"
  wc -l "$BAD"
done
