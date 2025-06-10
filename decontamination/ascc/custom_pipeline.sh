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

cd /scratch/leuven/357/vsc35707/decontamination/

# === Variables ===
THREADS=36
ASSEMBLY=purged.fa
KRAKEN_DB=./kraken2db
NR_DB=./nr/nr.dmnd
TAXDUMP=./taxdump
OUTPUT=decontam_output_regions
TARGET_TAXID=235516  # Pogonus chalceus

mkdir -p $OUTPUT

# === Step 1: Kraken2 classification ===
kraken2 --db $KRAKEN_DB --threads $THREADS \
  --report $OUTPUT/kraken2.report \
  --output $OUTPUT/kraken2.out \
  $ASSEMBLY

# Extract contigs assigned to non-target taxIDs (contaminants)
awk -v target=$TARGET_TAXID '$6 != target {print $2}' $OUTPUT/kraken2.out > $OUTPUT/kraken_contaminants.txt

# === Step 2: Tiara classification ===
tiara -i $ASSEMBLY -o $OUTPUT/tiara.tsv -t $THREADS
awk '$2 != "eukaryote"' $OUTPUT/tiara.tsv | cut -f1 > $OUTPUT/tiara_contaminants.txt

# === Step 3: DIAMOND BLASTX (nr) ===
diamond blastx -d $NR_DB -q $ASSEMBLY -o $OUTPUT/diamond.tsv \
  -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle \
  -k 1 --evalue 1e-5 --threads $THREADS

#wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
#gunzip prot.accession2taxid.gz
#cut -f1,3 prot.accession2taxid > acc2taxid.tsv

# Extract unique subject accessions from DIAMOND output
cut -f2 $OUTPUT/diamond.tsv | sort | uniq > $OUTPUT/diamond_subjects.txt

grep -Fwf $OUTPUT/diamond_subjects.txt acc2taxid.tsv > $OUTPUT/filtered_acc2taxid.tsv

# Extract only hits to contaminant taxids: 2 (bacteria), 4751 (fungi), 2157 (archaea), 9606 (human)
grep -Fwf <(awk '$2 ~ /^(2|4751|2157|9606)$/' $OUTPUT/filtered_acc2taxid.tsv | cut -f1) $OUTPUT/diamond.tsv > $OUTPUT/diamond_contaminants.tsv

# Convert DIAMOND contaminant hits to BED
awk 'BEGIN{OFS="\t"} {s=($7<$8)?$7:$8; e=($7>$8)?$7:$8; print $1, s-1, e}' $OUTPUT/diamond_contaminants.tsv > $OUTPUT/contaminant_regions.bed

# === Step 5: Extract contaminant sequences ===
bedtools getfasta -fi $ASSEMBLY -bed $OUTPUT/contaminant_regions.bed -fo $OUTPUT/contaminants_from_regions.fasta

# === Step 6: Mask contaminant regions ===
bedtools maskfasta -fi $ASSEMBLY -bed $OUTPUT/contaminant_regions.bed -fo $OUTPUT/assembly.masked.fasta

# === Step 7: Split at Ns and remove small fragments ===
# Split contigs at stretches of Ns (>10 Ns)
seqkit seq -w 0 $OUTPUT/assembly.masked.fasta > $OUTPUT/flat_masked.fasta
seqkit sliding -s 1 -W 10000 $OUTPUT/flat_masked.fasta | seqkit seq -m 500 > $OUTPUT/assembly.cleaned.fasta

# === Combine all contaminant IDs ===
cat $OUTPUT/kraken_contaminants.txt $OUTPUT/tiara_contaminants.txt <(cut -f1 $OUTPUT/diamond_contaminants.tsv) | sort | uniq > $OUTPUT/all_contaminant_contigs.txt

# === BUSCO with available lineages ===
LINEAGES=(archaea_odb12 bacteria_odb12 coleoptera_odb12 eukaryota_odb12 fungi_odb12 mammalia_odb12)
for lineage in "${LINEAGES[@]}"; do
  compleasm.py run -a $OUTPUT/assembly.cleaned.fasta -o $OUTPUT/busco_${lineage} -l $lineage -t $THREADS
done

echo "=== Decontamination complete. Cleaned assembly: $OUTPUT/assembly.cleaned.fasta ==="
