#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=final_fasta_all
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --time=4:00:00
#SBATCH -o final_fasta_all.%j.out
#SBATCH -A lp_edu_eeg_2026

LEPANCHOR="/data/leuven/357/vsc35707/LepAnchor"
GENOME="/scratch/leuven/357/vsc35707/sw-assembly/T2T_assembly/yahs_scaffolding/Pchalceus_SW_yahs_scaffolds_final.fa"

AGP="map15.agp"

OUT="SW.final.with_unplaced.fasta"
CHR_ONLY="SW.final.chromosomes_only.fasta"

# Make chromosomes from AGP
awk -f "$LEPANCHOR/makefasta.awk" "$GENOME" "$AGP" > "$CHR_ONLY"

# Reformat chromosomes + append unused scaffolds cleanly
python << EOF
from pathlib import Path
import re

genome = Path("$GENOME")
agp = Path("$AGP")
chr_only = Path("$CHR_ONLY")
out = Path("$OUT")

def read_fasta(path):
    seqs = {}
    order = []
    name = None
    chunks = []

    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).upper()
                    order.append(name)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)

    if name is not None:
        seqs[name] = "".join(chunks).upper()
        order.append(name)

    return seqs, order

def write_record(fout, name, seq, width=80):
    fout.write(f">{name}\\n")
    for i in range(0, len(seq), width):
        fout.write(seq[i:i+width] + "\\n")

# scaffolds used in AGP
used = set()
with agp.open() as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.rstrip().split()
        if len(fields) >= 9 and fields[4] == "W":
            used.add(fields[5])

chr_seqs, chr_order = read_fasta(chr_only)
genome_seqs, genome_order = read_fasta(genome)

def scaffold_number(x):
    m = re.search(r"scaffold_(\\d+)$", x)
    return int(m.group(1)) if m else 10**12

remaining = [s for s in genome_order if s not in used]
remaining = sorted(remaining, key=scaffold_number)

with out.open("w") as fout:
    for name in chr_order:
        write_record(fout, name, chr_seqs[name])

    for name in remaining:
        write_record(fout, name, genome_seqs[name])

print(f"Used scaffolds in AGP: {len(used)}")
print(f"Chromosomes written: {len(chr_order)}")
print(f"Unplaced scaffolds appended: {len(remaining)}")
print(f"Final FASTA: {out}")
EOF

echo "Done. Now index with a working samtools:"
echo "samtools faidx $OUT"
