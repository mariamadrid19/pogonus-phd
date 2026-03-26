#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=qc_rna_bams
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=08:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o qc_rna_bams.%j.out

# ------------------------------------------------------------
# QC for RNA-seq BAMs:
# - STAR Illumina BAMs
# - minimap2 Iso-Seq BAMs
#
# Outputs:
#   qc_bam_summary.tsv
#   qc_idxstats.tsv
# ------------------------------------------------------------

conda activate thesis 

# ---- User settings ----
STAR_DIR="/scratch/leuven/357/vsc35707/annotation/star_bams"
ISOSEQ_DIR="/scratch/leuven/357/vsc35707/annotation/isoseq_bams"
OUTDIR="/scratch/leuven/357/vsc35707/annotation/qc_rnaseq_bams"
THREADS=4

mkdir -p "$OUTDIR"

SUMMARY="${OUTDIR}/qc_bam_summary.tsv"
IDXSTATS="${OUTDIR}/qc_idxstats.tsv"

# ---- Header for summary table ----
printf "sample\ttype\tbam\tquickcheck\ttotal_reads\tmapped_reads\tmapped_pct\tsecondary\t%s\tduplicates\tpaired_reads\tproperly_paired\tproperly_paired_pct\tmapq0_reads\tmapq0_pct_
of_mapped\tspliced_reads\tspliced_pct_of_mapped\n" "supplementary" > "$SUMMARY"

# ---- Header for idxstats table ----
printf "sample\ttype\tchromosome\tchrom_length\tmapped_reads\tunmapped_reads\n" > "$IDXSTATS"

# ---- Collect BAMs ----
shopt -s nullglob
bams=(
  "${STAR_DIR}"/*.Aligned.sortedByCoord.out.bam
  "${ISOSEQ_DIR}"/*.sorted.bam
)

if [[ ${#bams[@]} -eq 0 ]]; then
  echo "ERROR: no BAM files found."
  exit 1
fi

for bam in "${bams[@]}"; do
  base=$(basename "$bam")

  # Infer sample and type
  if [[ "$base" == *.Aligned.sortedByCoord.out.bam ]]; then
    sample="${base%.Aligned.sortedByCoord.out.bam}"
    bam_type="Illumina"
  elif [[ "$base" == *.sorted.bam ]]; then
    sample="${base%.sorted.bam}"
    bam_type="IsoSeq"
  else
    sample="$base"
    bam_type="Unknown"
  fi

  echo "Processing: $sample ($bam_type)"

  # 1) Quick structural check
  if samtools quickcheck "$bam"; then
    quick="PASS"
  else
    quick="FAIL"
  fi

  # 2) flagstat summary
  flagfile="${OUTDIR}/${sample}.flagstat.txt"
  samtools flagstat "$bam" > "$flagfile"

  total_reads=$(awk '/in total/ {print $1; exit}' "$flagfile")
  secondary=$(awk '/secondary/ {print $1; exit}' "$flagfile")
  supplementary=$(awk '/supplementary/ {print $1; exit}' "$flagfile")
  duplicates=$(awk '/duplicates/ {print $1; exit}' "$flagfile")
  mapped_reads=$(awk '/mapped \(/ && $0 !~ /mate mapped/ {print $1; exit}' "$flagfile")
  paired_reads=$(awk '/paired in sequencing/ {print $1; exit}' "$flagfile")
  properly_paired=$(awk '/properly paired/ {print $1; exit}' "$flagfile")

  # 3) Percentages
  mapped_pct=$(awk -v m="$mapped_reads" -v t="$total_reads" 'BEGIN{if(t>0) printf "%.4f", 100*m/t; else print "NA"}')
  properly_paired_pct=$(awk -v p="$properly_paired" -v t="$paired_reads" 'BEGIN{if(t>0) printf "%.4f", 100*p/t; else print "NA"}')

  # 4) MAPQ 0 count among mapped alignments
  mapq0_reads=$(samtools view "$bam" | awk '$5==0{c++} END{print c+0}')
  mapq0_pct=$(awk -v q="$mapq0_reads" -v m="$mapped_reads" 'BEGIN{if(m>0) printf "%.4f", 100*q/m; else print "NA"}')

  # 5) Spliced alignments: CIGAR contains N
  spliced_reads=$(samtools view "$bam" | awk '$6 ~ /N/ {c++} END{print c+0}')
  spliced_pct=$(awk -v s="$spliced_reads" -v m="$mapped_reads" 'BEGIN{if(m>0) printf "%.4f", 100*s/m; else print "NA"}')

  # 6) Write summary row
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$sample" "$bam_type" "$bam" "$quick" "$total_reads" "$mapped_reads" "$mapped_pct" \
    "$secondary" "$supplementary" "$duplicates" "$paired_reads" "$properly_paired" \
    "$properly_paired_pct" "$mapq0_reads" "$mapq0_pct" "$spliced_reads" "$spliced_pct" \
    >> "$SUMMARY"

  # 7) idxstats by chromosome/scaffold
  samtools idxstats "$bam" | \
    awk -v s="$sample" -v t="$bam_type" 'BEGIN{OFS="\t"} {print s, t, $1, $2, $3, $4}' \
    >> "$IDXSTATS"
done

echo
echo "Done."
echo "Summary table : $SUMMARY"
echo "Idxstats table: $IDXSTATS"
