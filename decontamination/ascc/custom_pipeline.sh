#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --partition=batch_long
#SBATCH --job-name=decontamination
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=150G
#SBATCH --time=100:00:00
#SBATCH -o decon.%j.out
#SBATCH -A lp_svbelleghem

# === Activate conda + modules ===
conda activate decon_env
module load BLAST+/2.13.0-gompi-2022a
module load SAMtools/1.16.1-GCC-11.3.0
module load BEDTools/2.27.1-intel-2018a

cd /scratch/leuven/357/vsc35707/decontamination/

# === Variables ===
THREADS=36
ASSEMBLY=purged.fa
KRAKEN_DB=./kraken2db
NR_DB=./nr/nr.dmnd
OUTPUT=decontam_output_regions
TARGET_TAXID=235516  # Pogonus chalceus

mkdir -p $OUTPUT

# === Step 1: Kraken2 classification ===
kraken2 --db $KRAKEN_DB --threads $THREADS \
  --report $OUTPUT/kraken2.report \
  --output $OUTPUT/kraken2.out \
  $ASSEMBLY

# Extract non-target contigs (not Pogonus chalceus)
awk -v target=$TARGET_TAXID '$6 != target {print $2}' $OUTPUT/kraken2.out > $OUTPUT/kraken_contaminants.txt

# === Step 2: Tiara classification ===
tiara -i $ASSEMBLY -o $OUTPUT/tiara.tsv -t $THREADS
awk '$2 != "eukaryote"' $OUTPUT/tiara.tsv | cut -f1 > $OUTPUT/tiara_contaminants.txt

# === Step 3: DIAMOND BLASTX (nr) ===
diamond blastx -d $NR_DB -q $ASSEMBLY -o $OUTPUT/diamond.tsv \
  -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  -k 1 --evalue 1e-5 --threads $THREADS

# === Step 4: Map DIAMOND sseqids to taxids ===
cut -f2 $OUTPUT/diamond.tsv | sort | uniq > $OUTPUT/diamond_subjects.txt
grep -Fwf $OUTPUT/diamond_subjects.txt acc2taxid.tsv > $OUTPUT/filtered_acc2taxid.tsv

# Extract DIAMOND hits to contaminant taxids
grep -Fwf <(awk '$2 ~ /^(2|4751|2157|9606)$/' $OUTPUT/filtered_acc2taxid.tsv | cut -f1) $OUTPUT/diamond.tsv > $OUTPUT/diamond_contaminants.tsv

# Convert DIAMOND hits to BED
awk 'BEGIN{OFS="\t"} {s=($7<$8)?$7:$8; e=($7>$8)?$7:$8; print $1, s-1, e}' $OUTPUT/diamond_contaminants.tsv > $OUTPUT/contaminant_regions_diamond.bed

# === Step 5: Get full contig lengths ===
seqkit fx2tab -n -l $ASSEMBLY > $OUTPUT/contig_lengths.tsv

# Filter contig-lengths for Kraken2 + Tiara contaminant contigs
grep -Fwf $OUTPUT/kraken_contaminants.txt $OUTPUT/contig_lengths.tsv > $OUTPUT/kraken_contaminant_locs.tsv
grep -Fwf $OUTPUT/tiara_contaminants.txt $OUTPUT/contig_lengths.tsv > $OUTPUT/tiara_contaminant_locs.tsv

# Convert Kraken2 and Tiara hits to full BED spans
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' $OUTPUT/kraken_contaminant_locs.tsv > $OUTPUT/contaminant_regions_kraken.bed
awk 'BEGIN{OFS="\t"} {print $1, 0, $2}' $OUTPUT/tiara_contaminant_locs.tsv > $OUTPUT/contaminant_regions_tiara.bed

# === Step 6: Merge all regions into one BED ===
cat $OUTPUT/contaminant_regions_diamond.bed \
    $OUTPUT/contaminant_regions_kraken.bed \
    $OUTPUT/contaminant_regions_tiara.bed \
    | sort -k1,1 -k2,2n > $OUTPUT/all_contaminant_regions.bed

# Merge overlapping or adjacent regions
bedtools merge -i $OUTPUT/all_contaminant_regions.bed > $OUTPUT/merged_contaminant_regions.bed

# === Step 7: Extract contaminant sequences (for record) ===
bedtools getfasta -fi $ASSEMBLY -bed $OUTPUT/merged_contaminant_regions.bed -fo $OUTPUT/contaminants_from_regions.fasta

# === Step 8: Mask contaminant regions ===
bedtools maskfasta -fi $ASSEMBLY -bed $OUTPUT/merged_contaminant_regions.bed -fo $OUTPUT/assembly.masked.fasta

# === Step 9: Split masked contigs and remove tiny fragments ===
seqkit seq -w 0 $OUTPUT/assembly.masked.fasta > $OUTPUT/flat_masked.fasta
seqkit sliding -s 1 -W 10000 $OUTPUT/flat_masked.fasta | seqkit seq -m 500 > $OUTPUT/purged.decont.fasta

# === Step 10: Combine all contaminant contig names ===
cat $OUTPUT/kraken_contaminants.txt $OUTPUT/tiara_contaminants.txt <(cut -f1 $OUTPUT/diamond_contaminants.tsv) | sort | uniq > $OUTPUT/all_contaminant_contigs.txt

# === Step 11: BUSCO with available lineages ===
LINEAGES=(archaea_odb12 bacteria_odb12 coleoptera_odb12 eukaryota_odb12 fungi_odb12 mammalia_odb12)
for lineage in "${LINEAGES[@]}"; do
  compleasm.py run -a $OUTPUT/purged.decont.fasta -o $OUTPUT/busco_${lineage} -l $lineage -t $THREADS
done

echo "=== Decontamination complete. Cleaned assembly: $OUTPUT/assembly.cleaned.fasta ==="
