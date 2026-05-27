#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=lepanchor_filter_grid
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=4:00:00
#SBATCH -o lepanchor_filter_grid.%j.out
#SBATCH -A lp_edu_eeg_2026

module load Java/25.36
conda activate variant_tools

set -euo pipefail

cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap/results-map12-t05

LEPANCHOR="/data/leuven/357/vsc35707/LepAnchor"

MAP_PREFIX="map12"
NCHR=11
MINPUT_SUFFIX="_m"
MAP_BED="map.bed"

# =========================================================
# Run one filtering + LepAnchor analysis
# =========================================================

run_filter_set () {
  LABEL="$1"

  MIN_MARKERS_VERTICAL="$2"
  MIN_SPAN_BP_VERTICAL="$3"
  MIN_SPAN_FRAC_VERTICAL="$4"

  BP_BIN_SIZE="$5"
  MIN_MARKERS_HORIZONTAL="$6"
  MIN_CM_SPAN_HORIZONTAL="$7"

  OUTDIR="filter_${LABEL}"

  echo ""
  echo "===================================================="
  echo "Running filter set: ${LABEL}"
  echo "Output directory: ${OUTDIR}"
  echo "Vertical:"
  echo "  MIN_MARKERS_VERTICAL=${MIN_MARKERS_VERTICAL}"
  echo "  MIN_SPAN_BP_VERTICAL=${MIN_SPAN_BP_VERTICAL}"
  echo "  MIN_SPAN_FRAC_VERTICAL=${MIN_SPAN_FRAC_VERTICAL}"
  echo "Horizontal:"
  echo "  BP_BIN_SIZE=${BP_BIN_SIZE}"
  echo "  MIN_MARKERS_HORIZONTAL=${MIN_MARKERS_HORIZONTAL}"
  echo "  MIN_CM_SPAN_HORIZONTAL=${MIN_CM_SPAN_HORIZONTAL}"
  echo "===================================================="
  echo ""

  rm -rf "$OUTDIR"

  mkdir -p \
    "${OUTDIR}" \
    "${OUTDIR}/clean_reports" \
    "${OUTDIR}/clean_las" \
    "${OUTDIR}/clean_agp"

  python << EOF
import pandas as pd
import os

MAP_PREFIX = "${MAP_PREFIX}"
NCHR = ${NCHR}
MINPUT_SUFFIX = "${MINPUT_SUFFIX}"
OUTDIR = "${OUTDIR}"

MIN_MARKERS_VERTICAL = ${MIN_MARKERS_VERTICAL}
MIN_SPAN_BP_VERTICAL = ${MIN_SPAN_BP_VERTICAL}
MIN_SPAN_FRAC_VERTICAL = ${MIN_SPAN_FRAC_VERTICAL}

BP_BIN_SIZE = ${BP_BIN_SIZE}
MIN_MARKERS_HORIZONTAL = ${MIN_MARKERS_HORIZONTAL}
MIN_CM_SPAN_HORIZONTAL = ${MIN_CM_SPAN_HORIZONTAL}

report_dir = f"{OUTDIR}/clean_reports"

summary_rows = []

for chrom in range(1, NCHR + 1):
    infile = f"{MAP_PREFIX}_chr{chrom}{MINPUT_SUFFIX}.input"
    outfile = f"{OUTDIR}/{MAP_PREFIX}_chr{chrom}{MINPUT_SUFFIX}.clean.input"

    print(f"Cleaning chromosome {chrom}: {infile}")

    df = pd.read_csv(
        infile,
        sep=r"\\s+",
        header=None,
        names=["scaffold", "pos", "LG", "male_cM", "female_cM"],
        dtype={"scaffold": str}
    )

    for col in ["pos", "LG", "male_cM", "female_cM"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    malformed = df[df[["pos", "LG", "male_cM", "female_cM"]].isna().any(axis=1)]

    if len(malformed) > 0:
        malformed.to_csv(
            f"{report_dir}/{MAP_PREFIX}_chr{chrom}.malformed_rows.tsv",
            sep="\\t",
            index=False
        )
        print(f"  malformed rows dropped: {len(malformed)}")

    df = df.dropna(subset=["pos", "LG", "male_cM", "female_cM"]).copy()
    df["pos"] = df["pos"].astype(int)
    df["LG"] = df["LG"].astype(int)
    df["row_id"] = range(len(df))

    original_n = len(df)

    bad_vertical = set()
    bad_horizontal = set()

    # =====================================================
    # Vertical filter:
    # same cM, very different physical positions
    # =====================================================

    scaffold_span = df.groupby("scaffold")["pos"].agg(lambda x: x.max() - x.min())

    vertical_report_rows = []

    for key, g in df.groupby(["scaffold", "LG", "female_cM"]):
        scaffold, lg, female_cm = key

        span_bp = g["pos"].max() - g["pos"].min()
        total_span = scaffold_span.loc[scaffold]
        frac = span_bp / total_span if total_span > 0 else 0

        if (
            len(g) >= MIN_MARKERS_VERTICAL
            and span_bp >= MIN_SPAN_BP_VERTICAL
            and frac >= MIN_SPAN_FRAC_VERTICAL
        ):
            bad_vertical.update(g["row_id"].tolist())

            vertical_report_rows.append({
                "chromosome": chrom,
                "scaffold": scaffold,
                "LG": lg,
                "female_cM": female_cm,
                "n_markers": len(g),
                "min_bp": g["pos"].min(),
                "max_bp": g["pos"].max(),
                "span_bp": span_bp,
                "scaffold_marker_span": total_span,
                "span_fraction": frac
            })

    pd.DataFrame(vertical_report_rows).to_csv(
        f"{report_dir}/{MAP_PREFIX}_chr{chrom}.vertical_groups_report.tsv",
        sep="\\t",
        index=False
    )

    # =====================================================
    # Horizontal filter:
    # nearby physical position, very different cM
    # =====================================================

    df["bp_bin"] = (df["pos"] // BP_BIN_SIZE).astype(int)

    horizontal_report_rows = []

    for key, g in df.groupby(["scaffold", "LG", "bp_bin"]):
        scaffold, lg, bp_bin = key

        cm_span = g["female_cM"].max() - g["female_cM"].min()

        if (
            len(g) >= MIN_MARKERS_HORIZONTAL
            and cm_span >= MIN_CM_SPAN_HORIZONTAL
        ):
            bad_horizontal.update(g["row_id"].tolist())

            horizontal_report_rows.append({
                "chromosome": chrom,
                "scaffold": scaffold,
                "LG": lg,
                "bp_bin": bp_bin,
                "n_markers": len(g),
                "min_bp": g["pos"].min(),
                "max_bp": g["pos"].max(),
                "min_cM": g["female_cM"].min(),
                "max_cM": g["female_cM"].max(),
                "cm_span": cm_span
            })

    pd.DataFrame(horizontal_report_rows).to_csv(
        f"{report_dir}/{MAP_PREFIX}_chr{chrom}.horizontal_groups_report.tsv",
        sep="\\t",
        index=False
    )

    bad_all = bad_vertical | bad_horizontal

    keep_cols = ["scaffold", "pos", "LG", "male_cM", "female_cM"]

    df[df["row_id"].isin(bad_vertical)][keep_cols].to_csv(
        f"{report_dir}/{MAP_PREFIX}_chr{chrom}.bad_vertical_markers.tsv",
        sep="\\t",
        index=False
    )

    df[df["row_id"].isin(bad_horizontal)][keep_cols].to_csv(
        f"{report_dir}/{MAP_PREFIX}_chr{chrom}.bad_horizontal_markers.tsv",
        sep="\\t",
        index=False
    )

    df[df["row_id"].isin(bad_all)][keep_cols].to_csv(
        f"{report_dir}/{MAP_PREFIX}_chr{chrom}.bad_all_markers.tsv",
        sep="\\t",
        index=False
    )

    clean = df[~df["row_id"].isin(bad_all)][keep_cols]

    clean.to_csv(outfile, sep="\\t", header=False, index=False)

    summary_rows.append({
        "chromosome": chrom,
        "original_markers": original_n,
        "vertical_removed": len(bad_vertical),
        "horizontal_removed": len(bad_horizontal),
        "total_removed": len(bad_all),
        "cleaned_markers": len(clean),
        "percent_removed": round(100 * len(bad_all) / original_n, 2) if original_n > 0 else 0
    })

    print(f"  vertical removed:   {len(bad_vertical)}")
    print(f"  horizontal removed: {len(bad_horizontal)}")
    print(f"  total removed:      {len(bad_all)}")
    print(f"  remaining markers:  {len(clean)}")

pd.DataFrame(summary_rows).to_csv(
    f"{report_dir}/filter_summary.tsv",
    sep="\\t",
    index=False
)
EOF

  # =====================================================
  # Run PlaceAndOrientContigs
  # =====================================================

  for X in $(seq 1 "$NCHR"); do
    INPUT="${OUTDIR}/${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.clean.input"

    echo "Running PlaceAndOrientContigs for ${LABEL}, chromosome ${X}"

    java -cp "$LEPANCHOR/bin/" PlaceAndOrientContigs \
      map="$INPUT" \
      bed="$MAP_BED" \
      chromosome="$X" \
      noIntervals=1 \
      > "${OUTDIR}/clean_las/${MAP_PREFIX}_chr${X}.clean.la" \
      2> "${OUTDIR}/clean_las/${MAP_PREFIX}_chr${X}.clean.la.err"
  done

  # =====================================================
  # Generate AGP files
  # =====================================================

  for X in $(seq 1 "$NCHR"); do
    echo "Generating AGP for ${LABEL}, chromosome ${X}"

    awk -vlg="$X" \
      -f "$LEPANCHOR/makeagp_full2.awk" \
      "${OUTDIR}/clean_las/${MAP_PREFIX}_chr${X}.clean.la" \
      > "${OUTDIR}/clean_agp/${MAP_PREFIX}_chr${X}.clean.agp"
  done

  echo "Finished ${LABEL}. Results in ${OUTDIR}/"
}

# =========================================================
# Three filter strengths
# =========================================================

# LAX: remove only very obvious artifacts
run_filter_set \
  "lax" \
  25 15000000 0.50 \
  50000 10 40

# MEDIUM: balanced
run_filter_set \
  "medium" \
  15 10000000 0.35 \
  50000 5 30

# STRICT: remove more visible line artifacts
run_filter_set \
  "strict" \
  10 5000000 0.25 \
  25000 3 20

echo ""
echo "All runs complete."
echo "Plot these folders separately in R:"
echo "  filter_lax/"
echo "  filter_medium/"
echo "  filter_strict/"
