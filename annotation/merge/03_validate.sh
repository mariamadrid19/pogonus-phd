#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=validate_LW_merge
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH -A lp_svbelleghem
#SBATCH -o validate_LW_merge.%j.out
#SBATCH -e validate_LW_merge.%j.err

set -euo pipefail

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate agat-merge

output_dir=/scratch/leuven/357/vsc35707/annotation/LW-annotation/merge-annotations

merged="${output_dir}/Pchalceus_LW.BRAKER_EviAnn.merged.gff3"
braker="${output_dir}/braker.normalized.gff3"
eviann="${output_dir}/eviann.normalized.gff3"

# CHANGE THIS if the assembly is located elsewhere.
genome=/scratch/leuven/357/vsc35707/annotation/LW-annotation/Pchalceus_LW_final.fasta

test -s "$merged"
test -s "$braker"
test -s "$eviann"
test -s "$genome"

cd "$output_dir"

# 1. Confirm that all annotation scaffold names exist in the FASTA.
#samtools faidx "$genome"

cut -f1 "${genome}.fai" |
sort -u > genome.contigs

awk -F '\t' '$0 !~ /^#/ {print $1}' "$merged" |
sort -u > merged.contigs

comm -23 merged.contigs genome.contigs \
    > annotation_contigs_missing_from_fasta.txt

if [[ -s annotation_contigs_missing_from_fasta.txt ]]; then
    echo "ERROR: annotation contains sequence names absent from FASTA"
    cat annotation_contigs_missing_from_fasta.txt
    exit 1
fi

echo "All annotation sequence names occur in the FASTA"

# 2. Check that every Parent identifier exists.
awk -F '\t' '
$0 !~ /^#/ && NF == 9 {
    n=split($9,attributes,";")

    for (i=1; i<=n; i++) {
        if (attributes[i] ~ /^ID=/) {
            id=attributes[i]
            sub(/^ID=/,"",id)
            ids[id]=1
        }

        if (attributes[i] ~ /^Parent=/) {
            parent_value=attributes[i]
            sub(/^Parent=/,"",parent_value)
            number_parents=split(parent_value,parent_parts,",")

            for (j=1; j<=number_parents; j++)
                parents[parent_parts[j]]=1
        }
    }
}
END {
    for (parent_id in parents)
        if (!(parent_id in ids))
            print parent_id
}' "$merged" |
sort -u > missing_parent_ids.txt



if [[ -s missing_parent_ids.txt ]]; then
    echo "ERROR: missing parent IDs found"
    head -20 missing_parent_ids.txt
    exit 1
fi

echo "All Parent identifiers resolve"

# 3. Extract merged protein sequences.
# --clean_final_stop removes only the normal terminal '*'.
# Internal stop codons remain visible.
agat_sp_extract_sequences.pl \
    --gff "$merged" \
    --fasta "$genome" \
    --type cds \
    --protein \
    --clean_final_stop \
    --output Pchalceus_LW.merged.proteins.faa

# 4. Extract CDS nucleotide sequences.
agat_sp_extract_sequences.pl \
    --gff "$merged" \
    --fasta "$genome" \
    --type cds \
    --output Pchalceus_LW.merged.cds.fna

# 5. Report proteins containing internal stop codons.
awk '
function finish_sequence() {
    if (header != "" && sequence ~ /\*/)
        print header
}
/^>/ {
    finish_sequence()
    header=$0
    sequence=""
    next
}
{
    sequence=sequence $0
}
END {
    finish_sequence()
}' Pchalceus_LW.merged.proteins.faa \
> proteins_with_internal_stops.txt

# 6. Report protein lengths.
awk '
function finish_sequence() {
    if (header != "")
        print header "\t" length(sequence)
}
/^>/ {
    finish_sequence()
    header=substr($0,2)
    sequence=""
    next
}
{
    sequence=sequence $0
}
END {
    finish_sequence()
}' Pchalceus_LW.merged.proteins.faa \
> merged.protein_lengths.tsv

awk -F '\t' '$2 < 50' merged.protein_lengths.tsv \
    > proteins_shorter_than_50aa.tsv

awk -F '\t' '$2 > 10000' merged.protein_lengths.tsv \
    > proteins_longer_than_10000aa.tsv

echo "Protein sequences:"
grep -c '^>' Pchalceus_LW.merged.proteins.faa

echo "Proteins with internal stop codons:"
wc -l < proteins_with_internal_stops.txt

echo "Proteins shorter than 50 aa:"
wc -l < proteins_shorter_than_50aa.tsv

echo "Proteins longer than 10,000 aa:"
wc -l < proteins_longer_than_10000aa.tsv

# 7. Record attribute types, including functional annotations.
agat_sq_list_attributes.pl \
    --gff "$braker" \
    --output braker.attribute_inventory.txt

agat_sq_list_attributes.pl \
    --gff "$merged" \
    --output merged.attribute_inventory.txt

echo "Validation completed"
