#!/bin/bash
#SBATCH --cluster=wice
#SBATCH --job-name=LW_SW_gene_matching
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=LW_SW_gene_matching.%j.out
#SBATCH --error=LW_SW_gene_matching.%j.err
#SBATCH -A lp_svbelleghem

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate mmseqs2

LW_AA="/scratch/leuven/357/vsc35707/annotation/LW-annotation/functional/results/01_longest_isoforms/braker.longest_isoforms.aa"
SW_AA="/scratch/leuven/357/vsc35707/annotation/SW-annotation/functional/results/01_longest_isoforms/braker.longest_isoforms.aa"

OUTDIR="/scratch/leuven/357/vsc35707/rna-seq/LW_SW_gene_comparison"

mkdir -p "${OUTDIR}"
mkdir -p "${OUTDIR}/tmp_rbh"
mkdir -p "${OUTDIR}/tmp_all"

# Convert headers such as g1.t1 to gene IDs such as g1.
# Also remove terminal * characters from protein sequences.

awk '
/^>/ {
    header = $0
    sub(/^>[[:space:]]*/, "", header)
    split(header, fields, /[[:space:]]+/)
    id = fields[1]
    sub(/\.t[0-9]+$/, "", id)
    print ">" id
    next
}
{
    gsub(/\*/, "")
    print
}
' "${LW_AA}" > "${OUTDIR}/LW_longest_by_gene.aa"

awk '
/^>/ {
    header = $0
    sub(/^>[[:space:]]*/, "", header)
    split(header, fields, /[[:space:]]+/)
    id = fields[1]
    sub(/\.t[0-9]+$/, "", id)
    print ">" id
    next
}
{
    gsub(/\*/, "")
    print
}
' "${SW_AA}" > "${OUTDIR}/SW_longest_by_gene.aa"

# Verify that there are no duplicated gene IDs after removing .t1, .t2, etc.

grep '^>' "${OUTDIR}/LW_longest_by_gene.aa" |
    sed 's/^>//' |
    sort |
    uniq -d > "${OUTDIR}/LW_duplicate_gene_ids.txt"

grep '^>' "${OUTDIR}/SW_longest_by_gene.aa" |
    sed 's/^>//' |
    sort |
    uniq -d > "${OUTDIR}/SW_duplicate_gene_ids.txt"

echo "LW duplicated IDs:"
wc -l "${OUTDIR}/LW_duplicate_gene_ids.txt"

echo "SW duplicated IDs:"
wc -l "${OUTDIR}/SW_duplicate_gene_ids.txt"

# Reciprocal best-hit comparison

mmseqs easy-rbh \
    "${OUTDIR}/LW_longest_by_gene.aa" \
    "${OUTDIR}/SW_longest_by_gene.aa" \
    "${OUTDIR}/LW_vs_SW_rbh.raw.tsv" \
    "${OUTDIR}/tmp_rbh" \
    --threads 8 \
    --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits"

# Add descriptive column names

{
    printf "LW_gene_id\tSW_gene_id\tpident\talnlen\tqcov\ttcov\tevalue\tbits\n"
    cat "${OUTDIR}/LW_vs_SW_rbh.raw.tsv"
} > "${OUTDIR}/LW_vs_SW_rbh.tsv"

# Retain high-confidence one-to-one matches:
# at least 90% amino-acid identity and 80% coverage of both proteins.

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
NR == 1 {
    print
    next
}
$3 >= 90 && $5 >= 0.80 && $6 >= 0.80 {
    print
}
' "${OUTDIR}/LW_vs_SW_rbh.tsv" \
> "${OUTDIR}/LW_vs_SW_rbh.high_confidence.tsv"

# All-versus-all search for possible split, merged or duplicated genes

mmseqs easy-search \
    "${OUTDIR}/LW_longest_by_gene.aa" \
    "${OUTDIR}/SW_longest_by_gene.aa" \
    "${OUTDIR}/LW_vs_SW_all.raw.tsv" \
    "${OUTDIR}/tmp_all" \
    --threads 8 \
    -e 1e-10 \
    --format-output "query,target,pident,alnlen,qcov,tcov,evalue,bits"

{
    printf "LW_gene_id\tSW_gene_id\tpident\talnlen\tqcov\ttcov\tevalue\tbits\n"
    cat "${OUTDIR}/LW_vs_SW_all.raw.tsv"
} > "${OUTDIR}/LW_vs_SW_all.tsv"

# More permissive candidates for investigating non-one-to-one matches

awk -F '\t' '
BEGIN {
    OFS = "\t"
}
NR == 1 {
    print
    next
}
$3 >= 80 && $5 >= 0.70 && $6 >= 0.70 {
    print
}
' "${OUTDIR}/LW_vs_SW_all.tsv" \
> "${OUTDIR}/LW_vs_SW_all.filtered.tsv"

echo "Number of LW proteins:"
grep -c '^>' "${OUTDIR}/LW_longest_by_gene.aa"

echo "Number of SW proteins:"
grep -c '^>' "${OUTDIR}/SW_longest_by_gene.aa"

echo "Number of reciprocal best hits:"
wc -l "${OUTDIR}/LW_vs_SW_rbh.raw.tsv"

echo "Number of high-confidence reciprocal best hits:"
tail -n +2 "${OUTDIR}/LW_vs_SW_rbh.high_confidence.tsv" | wc -l

echo "Results written to ${OUTDIR}"
