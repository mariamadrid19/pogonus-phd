#!/bin/bash
#SBATCH --cluster=wice
#SBATCH --job-name=LW_SW_position_validation
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=LW_SW_position_validation.%j.out
#SBATCH --error=LW_SW_position_validation.%j.err
#SBATCH -A lp_svbelleghem

set -euo pipefail

# ============================================================
# Load software
# ============================================================

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate variant_tools

# ============================================================
# Input files
# ============================================================

# LW genome is the reference
LW_GENOME="/scratch/leuven/357/vsc35707/annotation/LW-annotation/Pchalceus_LW_final.fasta"

# SW genome is aligned and lifted onto LW
SW_GENOME="/scratch/leuven/357/vsc35707/annotation/SW-annotation/Pchalceus_SW_LW_chromosome_names.fasta"

# BRAKER GTF files
LW_GTF="/scratch/leuven/357/vsc35707/annotation/LW-annotation/braker-LW/braker.gtf"

SW_GTF="/scratch/leuven/357/vsc35707/annotation/SW-annotation/braker-SW/braker.gtf"

# High-confidence MMseqs2 reciprocal-best-hit table
RBH_FILE="/scratch/leuven/357/vsc35707/rna-seq/LW_SW_gene_comparison/LW_vs_SW_rbh.high_confidence.tsv"

# Output directory
OUTDIR="/scratch/leuven/357/vsc35707/rna-seq/LW_SW_position_validation"

# ============================================================
# Parameters
# ============================================================

THREADS=16

# Maximum distance between non-overlapping genes for "nearby"
NEARBY_DISTANCE=10000

# Minimum fraction of either gene interval that must overlap
MIN_OVERLAP_FRACTION=0.50

# Minimum assembly-alignment block accepted by paftools liftover
LIFTOVER_MIN_ALIGNMENT_LENGTH=1000

# Minimum mapping quality accepted by paftools liftover
LIFTOVER_MIN_MAPQ=5

# Reuse existing minimap2 PAF
FORCE_REALIGN=0

# Rerun liftover using the corrected settings
FORCE_LIFTOVER=1

# ============================================================
# Create output directory
# ============================================================

mkdir -p "${OUTDIR}"

# ============================================================
# Output file names
# ============================================================

LW_BED="${OUTDIR}/LW_genes.bed"
SW_BED="${OUTDIR}/SW_genes.bed"

PAF="${OUTDIR}/SW_to_LW.asm5.paf"
SORTED_PAF="${OUTDIR}/SW_to_LW.asm5.sorted.paf"
PAF_STATS="${OUTDIR}/SW_to_LW.asm5.stats.txt"

LIFTED_SW_BED="${OUTDIR}/SW_genes_lifted_to_LW.bed"

LIFTED_SW_BED_WITH_IDS="${OUTDIR}/SW_genes_lifted_to_LW.with_gene_ids.bed"

UNMATCHED_LIFTOVER_IDS="${OUTDIR}/SW_lifted_coordinate_ids_unmatched.tsv"

VALIDATION_TABLE="${OUTDIR}/LW_SW_RBH_position_validation.tsv"

SUMMARY_TABLE="${OUTDIR}/LW_SW_RBH_position_summary.tsv"

SUPPORTED_TABLE="${OUTDIR}/LW_SW_RBH_protein_and_position_supported.tsv"

QUESTIONABLE_TABLE="${OUTDIR}/LW_SW_RBH_position_questionable.tsv"

# ============================================================
# Check input files
# ============================================================

for FILE in \
    "${LW_GENOME}" \
    "${SW_GENOME}" \
    "${LW_GTF}" \
    "${SW_GTF}" \
    "${RBH_FILE}"
do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Input file is missing or empty:"
        echo "${FILE}"
        exit 1
    fi
done

for PROGRAM in minimap2 paftools.js gawk
do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} is not available."
        exit 1
    fi
done

# ============================================================
# Report inputs
# ============================================================

echo "Starting LW/SW genomic-position validation"
echo "LW genome: ${LW_GENOME}"
echo "SW genome: ${SW_GENOME}"
echo "LW GTF: ${LW_GTF}"
echo "SW GTF: ${SW_GTF}"
echo "RBH table: ${RBH_FILE}"
echo "Output directory: ${OUTDIR}"
echo "Minimap2 version: $(minimap2 --version)"
echo "GNU awk version: $(gawk --version | head -1)"
echo "Liftover minimum alignment length: ${LIFTOVER_MIN_ALIGNMENT_LENGTH}"
echo "Liftover minimum mapping quality: ${LIFTOVER_MIN_MAPQ}"

# ============================================================
# 1. Extract one genomic interval per BRAKER gene
# ============================================================

extract_gene_bed() {
    local GTF_FILE="$1"
    local BED_FILE="$2"

    gawk -F '\t' '
    BEGIN {
        OFS = "\t"
    }

    /^#/ || NF < 9 {
        next
    }

    {
        gene = ""
        number_attributes = split($9, attributes, ";")

        for (i = 1; i <= number_attributes; i++) {
            attribute = attributes[i]
            sub(/^[[:space:]]+/, "", attribute)

            if (attribute ~ /^gene_id[[:space:]]/) {
                sub(/^gene_id[[:space:]]+/, "", attribute)
                gsub(/"/, "", attribute)
                sub(/[[:space:]].*$/, "", attribute)

                gene = attribute
                break
            }

            if (attribute ~ /^gene_id=/) {
                sub(/^gene_id=/, "", attribute)
                gsub(/"/, "", attribute)
                sub(/[[:space:]].*$/, "", attribute)

                gene = attribute
                break
            }
        }

        if (gene == "") {
            next
        }

        chromosome[gene] = $1
        strand[gene] = $7

        # Convert GTF to BED coordinates
        bed_start = $4 - 1
        bed_end = $5

        if (!(gene in minimum_start) || bed_start < minimum_start[gene]) {
            minimum_start[gene] = bed_start
        }

        if (!(gene in maximum_end) || bed_end > maximum_end[gene]) {
            maximum_end[gene] = bed_end
        }
    }

    END {
        for (gene in chromosome) {
            print chromosome[gene], \
                  minimum_start[gene], \
                  maximum_end[gene], \
                  gene, \
                  0, \
                  strand[gene]
        }
    }
    ' "${GTF_FILE}" |
        sort -k1,1 -k2,2n > "${BED_FILE}"
}

extract_gene_bed \
    "${LW_GTF}" \
    "${LW_BED}"

extract_gene_bed \
    "${SW_GTF}" \
    "${SW_BED}"

echo "Number of LW genes extracted:"
wc -l "${LW_BED}"

echo "Number of SW genes extracted:"
wc -l "${SW_BED}"

if [[ ! -s "${LW_BED}" ]]; then
    echo "ERROR: No LW gene coordinates were extracted."
    exit 1
fi

if [[ ! -s "${SW_BED}" ]]; then
    echo "ERROR: No SW gene coordinates were extracted."
    exit 1
fi

# ============================================================
# 2. Check GTF and FASTA sequence names
# ============================================================

gawk '
/^>/ {
    name = substr($1, 2)
    print name
}
' "${LW_GENOME}" |
    sort -u > "${OUTDIR}/LW_fasta_sequence_names.txt"

gawk '
/^>/ {
    name = substr($1, 2)
    print name
}
' "${SW_GENOME}" |
    sort -u > "${OUTDIR}/SW_fasta_sequence_names.txt"

cut -f1 "${LW_BED}" |
    sort -u > "${OUTDIR}/LW_gtf_sequence_names.txt"

cut -f1 "${SW_BED}" |
    sort -u > "${OUTDIR}/SW_gtf_sequence_names.txt"

comm -23 \
    "${OUTDIR}/LW_gtf_sequence_names.txt" \
    "${OUTDIR}/LW_fasta_sequence_names.txt" \
    > "${OUTDIR}/LW_gtf_sequences_missing_from_fasta.txt"

comm -23 \
    "${OUTDIR}/SW_gtf_sequence_names.txt" \
    "${OUTDIR}/SW_fasta_sequence_names.txt" \
    > "${OUTDIR}/SW_gtf_sequences_missing_from_fasta.txt"

LW_MISSING_SEQUENCES=$(
    wc -l < "${OUTDIR}/LW_gtf_sequences_missing_from_fasta.txt"
)

SW_MISSING_SEQUENCES=$(
    wc -l < "${OUTDIR}/SW_gtf_sequences_missing_from_fasta.txt"
)

echo "LW GTF sequence names missing from LW FASTA: ${LW_MISSING_SEQUENCES}"
echo "SW GTF sequence names missing from SW FASTA: ${SW_MISSING_SEQUENCES}"

if [[ "${LW_MISSING_SEQUENCES}" -gt 0 ]]; then
    echo "WARNING: LW GTF sequences absent from LW FASTA:"
    head "${OUTDIR}/LW_gtf_sequences_missing_from_fasta.txt"
fi

if [[ "${SW_MISSING_SEQUENCES}" -gt 0 ]]; then
    echo "WARNING: SW GTF sequences absent from SW FASTA:"
    head "${OUTDIR}/SW_gtf_sequences_missing_from_fasta.txt"
fi

# ============================================================
# 3. Align SW genome to LW genome
#
# Reference: LW
# Query: SW
# ============================================================

if [[ "${FORCE_REALIGN}" -eq 1 || ! -s "${PAF}" ]]; then
    echo "Running minimap2 assembly alignment."

    minimap2 \
        -cx asm5 \
        --cs \
        --secondary=no \
        -t "${THREADS}" \
        "${LW_GENOME}" \
        "${SW_GENOME}" \
        > "${PAF}"
else
    echo "Existing minimap2 PAF found."
    echo "Skipping genome alignment: ${PAF}"
fi

if [[ ! -s "${PAF}" ]]; then
    echo "ERROR: minimap2 PAF file is empty."
    exit 1
fi

echo "Number of assembly-alignment records:"
wc -l "${PAF}"

sort \
    -k6,6 \
    -k8,8n \
    "${PAF}" \
    > "${SORTED_PAF}"

paftools.js stat \
    "${SORTED_PAF}" \
    > "${PAF_STATS}"

# ============================================================
# 4. Lift SW gene coordinates to LW coordinates
#
# paftools defaults to a minimum alignment length of 50 kb.
# Here it is reduced to 1 kb to retain shorter assembly blocks.
# ============================================================

if [[ "${FORCE_LIFTOVER}" -eq 1 || ! -s "${LIFTED_SW_BED}" ]]; then
    echo "Running paftools liftover."

    paftools.js liftover \
        -q "${LIFTOVER_MIN_MAPQ}" \
        -l "${LIFTOVER_MIN_ALIGNMENT_LENGTH}" \
        "${SORTED_PAF}" \
        "${SW_BED}" \
        > "${LIFTED_SW_BED}"
else
    echo "Existing lifted SW BED found."
    echo "Skipping liftover: ${LIFTED_SW_BED}"
fi

if [[ ! -s "${LIFTED_SW_BED}" ]]; then
    echo "ERROR: No SW genes were lifted to LW coordinates."
    exit 1
fi

echo "Number of raw lifted SW intervals:"
wc -l "${LIFTED_SW_BED}"

# ============================================================
# 5. Restore SW gene IDs after paftools liftover
#
# paftools replaces BED column 4 with:
#
#   query_sequence_start_end
#
# The original SW BED coordinates are used to translate this
# coordinate identifier back to the BRAKER gene ID.
# ============================================================

rm -f "${UNMATCHED_LIFTOVER_IDS}"

gawk \
    -v unmatched_file="${UNMATCHED_LIFTOVER_IDS}" \
    -F '\t' '
BEGIN {
    OFS = "\t"
}

# Read original SW gene coordinates
ARGIND == 1 {
    coordinate_id = $1 "_" $2 "_" $3
    gene_id[coordinate_id] = $4
    next
}

# Read raw paftools output
ARGIND == 2 {
    coordinate_id = $4

    # Remove clipping suffixes added by paftools
    sub(/_t5_t3$/, "", coordinate_id)
    sub(/_t5$/, "", coordinate_id)
    sub(/_t3$/, "", coordinate_id)

    if (coordinate_id in gene_id) {
        $4 = gene_id[coordinate_id]
        print
        matched++
    } else {
        print $0 > unmatched_file
        unmatched++
    }
}

END {
    print "Lifted records restored to gene IDs: " matched + 0 > "/dev/stderr"
    print "Lifted records without matching IDs: " unmatched + 0 > "/dev/stderr"
}
' \
    "${SW_BED}" \
    "${LIFTED_SW_BED}" \
    > "${LIFTED_SW_BED_WITH_IDS}"

if [[ ! -s "${LIFTED_SW_BED_WITH_IDS}" ]]; then
    echo "ERROR: No lifted intervals could be restored to SW gene IDs."
    exit 1
fi

echo "Number of lifted records with restored SW gene IDs:"
wc -l "${LIFTED_SW_BED_WITH_IDS}"

LIFTED_UNIQUE_GENES=$(
    cut -f4 "${LIFTED_SW_BED_WITH_IDS}" |
        sort -u |
        wc -l
)

echo "Number of unique SW genes lifted: ${LIFTED_UNIQUE_GENES}"

# ============================================================
# 6. Check how many RBH genes were lifted
# ============================================================

tail -n +2 "${RBH_FILE}" |
    cut -f2 |
    sort -u > "${OUTDIR}/RBH_SW_gene_ids.txt"

cut -f4 "${LIFTED_SW_BED_WITH_IDS}" |
    sort -u > "${OUTDIR}/lifted_SW_gene_ids.txt"

RBH_SW_COUNT=$(
    wc -l < "${OUTDIR}/RBH_SW_gene_ids.txt"
)

RBH_SW_LIFTED_COUNT=$(
    comm -12 \
        "${OUTDIR}/RBH_SW_gene_ids.txt" \
        "${OUTDIR}/lifted_SW_gene_ids.txt" |
        wc -l
)

echo "High-confidence RBH SW genes: ${RBH_SW_COUNT}"
echo "High-confidence RBH SW genes lifted: ${RBH_SW_LIFTED_COUNT}"

# ============================================================
# 7. Validate MMseqs2 matches using genomic coordinates
# ============================================================

gawk \
    -v nearby_distance="${NEARBY_DISTANCE}" \
    -v minimum_overlap_fraction="${MIN_OVERLAP_FRACTION}" \
    -F '\t' '
BEGIN {
    OFS = "\t"
}

# Read LW gene coordinates
ARGIND == 1 {
    lw_chr[$4] = $1
    lw_start[$4] = $2
    lw_end[$4] = $3
    lw_strand[$4] = $6
    next
}

# Read lifted SW coordinates with restored gene IDs
ARGIND == 2 {
    sw_gene = $4

    lifted_count[sw_gene]++
    key = sw_gene SUBSEP lifted_count[sw_gene]

    lifted_chr[key] = $1
    lifted_start[key] = $2
    lifted_end[key] = $3
    lifted_strand[key] = (NF >= 6 ? $6 : ".")

    next
}

# Print header
ARGIND == 3 && FNR == 1 {
    print "LW_gene_id", \
          "SW_gene_id", \
          "pident", \
          "qcov", \
          "tcov", \
          "LW_chromosome", \
          "LW_start", \
          "LW_end", \
          "SW_lifted_chromosome", \
          "SW_lifted_start", \
          "SW_lifted_end", \
          "overlap_bp", \
          "overlap_fraction_LW", \
          "overlap_fraction_SW", \
          "distance_bp", \
          "position_status", \
          "number_of_SW_liftovers"

    next
}

# Process MMseqs2 reciprocal-best-hit table
ARGIND == 3 {
    lw_gene = $1
    sw_gene = $2

    pident = $3
    qcov = $5
    tcov = $6

    if (!(lw_gene in lw_chr)) {
        print lw_gene, \
              sw_gene, \
              pident, \
              qcov, \
              tcov, \
              "NA", \
              "NA", \
              "NA", \
              "NA", \
              "NA", \
              "NA", \
              "NA", \
              "NA", \
              "NA", \
              "NA", \
              "LW_gene_missing_from_GTF", \
              lifted_count[sw_gene] + 0

        next
    }

    if (!(sw_gene in lifted_count)) {
        print lw_gene, \
              sw_gene, \
              pident, \
              qcov, \
              tcov, \
              lw_chr[lw_gene], \
              lw_start[lw_gene], \
              lw_end[lw_gene], \
              "NA", \
              "NA", \
              "NA", \
              0, \
              0, \
              0, \
              "NA", \
              "SW_not_lifted", \
              0

        next
    }

    for (i = 1; i <= lifted_count[sw_gene]; i++) {
        key = sw_gene SUBSEP i

        overlap = 0
        distance = "NA"
        status = ""

        lw_length = lw_end[lw_gene] - lw_start[lw_gene]
        sw_length = lifted_end[key] - lifted_start[key]

        overlap_fraction_lw = 0
        overlap_fraction_sw = 0

        if (lifted_chr[key] != lw_chr[lw_gene]) {
            status = "different_chromosome"
        } else {
            overlap_start = (lw_start[lw_gene] > lifted_start[key] ? lw_start[lw_gene] : lifted_start[key])
            overlap_end = (lw_end[lw_gene] < lifted_end[key] ? lw_end[lw_gene] : lifted_end[key])

            overlap = overlap_end - overlap_start

            if (overlap < 0) {
                overlap = 0
            }

            if (overlap > 0) {
                distance = 0
            } else if (lifted_end[key] <= lw_start[lw_gene]) {
                distance = lw_start[lw_gene] - lifted_end[key]
            } else {
                distance = lifted_start[key] - lw_end[lw_gene]
            }

            overlap_fraction_lw = (lw_length > 0 ? overlap / lw_length : 0)
            overlap_fraction_sw = (sw_length > 0 ? overlap / sw_length : 0)

            if (overlap_fraction_lw >= minimum_overlap_fraction || overlap_fraction_sw >= minimum_overlap_fraction) {
                status = "position_supported"
            } else if (overlap > 0) {
                status = "position_overlap_weak"
            } else if (distance <= nearby_distance) {
                status = "position_nearby"
            } else {
                status = "position_conflict"
            }
        }

        print lw_gene, \
              sw_gene, \
              pident, \
              qcov, \
              tcov, \
              lw_chr[lw_gene], \
              lw_start[lw_gene], \
              lw_end[lw_gene], \
              lifted_chr[key], \
              lifted_start[key], \
              lifted_end[key], \
              overlap, \
              overlap_fraction_lw, \
              overlap_fraction_sw, \
              distance, \
              status, \
              lifted_count[sw_gene]
    }
}
' \
    "${LW_BED}" \
    "${LIFTED_SW_BED_WITH_IDS}" \
    "${RBH_FILE}" \
    > "${VALIDATION_TABLE}"

if [[ ! -s "${VALIDATION_TABLE}" ]]; then
    echo "ERROR: Position-validation table is empty."
    exit 1
fi

# ============================================================
# 8. Produce position-status summary
# ============================================================

gawk -F '\t' '
BEGIN {
    OFS = "\t"
}

NR > 1 {
    count[$16]++
}

END {
    print "position_status", "number_of_records"

    for (status in count) {
        print status, count[status]
    }
}
' "${VALIDATION_TABLE}" > "${SUMMARY_TABLE}.unsorted"

{
    head -1 "${SUMMARY_TABLE}.unsorted"
    tail -n +2 "${SUMMARY_TABLE}.unsorted" |
        sort -k2,2nr
} > "${SUMMARY_TABLE}"

rm -f "${SUMMARY_TABLE}.unsorted"

# ============================================================
# 9. Extract protein and position-supported matches
# ============================================================

gawk -F '\t' '
BEGIN {
    OFS = "\t"
}

NR == 1 {
    print
    next
}

$16 == "position_supported" {
    print
}
' "${VALIDATION_TABLE}" > "${SUPPORTED_TABLE}"

# ============================================================
# 10. Extract questionable matches
# ============================================================

gawk -F '\t' '
BEGIN {
    OFS = "\t"
}

NR == 1 {
    print
    next
}

$16 != "position_supported" {
    print
}
' "${VALIDATION_TABLE}" > "${QUESTIONABLE_TABLE}"

# ============================================================
# 11. Approximate assembly-alignment statistics
# ============================================================

gawk -F '\t' '
BEGIN {
    OFS = "\t"
}

{
    query_length[$1] = $2
    aligned_blocks += $11
}

END {
    represented_query_length = 0

    for (sequence in query_length) {
        represented_query_length += query_length[sequence]
    }

    print "metric", "value"
    print "sum_PAF_alignment_block_lengths", aligned_blocks
    print "length_of_represented_SW_sequences", represented_query_length

    if (represented_query_length > 0) {
        print "approximate_alignment_ratio", aligned_blocks / represented_query_length
    } else {
        print "approximate_alignment_ratio", "NA"
    }
}
' "${PAF}" > "${OUTDIR}/SW_to_LW.approximate_alignment_summary.tsv"

# ============================================================
# 12. Final report
# ============================================================

SUPPORTED_COUNT=$(
    tail -n +2 "${SUPPORTED_TABLE}" |
        wc -l
)

QUESTIONABLE_COUNT=$(
    tail -n +2 "${QUESTIONABLE_TABLE}" |
        wc -l
)

echo
echo "Position-validation summary:"

if command -v column >/dev/null 2>&1; then
    column -t "${SUMMARY_TABLE}"
else
    cat "${SUMMARY_TABLE}"
fi

echo
echo "Input SW genes: $(wc -l < "${SW_BED}")"
echo "Unique SW genes lifted: ${LIFTED_UNIQUE_GENES}"
echo "High-confidence RBH SW genes: ${RBH_SW_COUNT}"
echo "High-confidence RBH SW genes lifted: ${RBH_SW_LIFTED_COUNT}"
echo "Position-supported RBH records: ${SUPPORTED_COUNT}"
echo "Questionable RBH records: ${QUESTIONABLE_COUNT}"

echo
echo "Main validation table:"
echo "${VALIDATION_TABLE}"

echo
echo "Protein and position-supported matches:"
echo "${SUPPORTED_TABLE}"

echo
echo "Questionable matches:"
echo "${QUESTIONABLE_TABLE}"

echo
echo "Assembly-alignment statistics:"
echo "${PAF_STATS}"

echo
echo "Approximate alignment summary:"
echo "${OUTDIR}/SW_to_LW.approximate_alignment_summary.tsv"

echo
echo "Finished successfully."
