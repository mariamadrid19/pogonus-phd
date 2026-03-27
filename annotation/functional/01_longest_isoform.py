#!/usr/bin/env python3

import argparse
import re
import sys
from collections import defaultdict


def parse_fasta(path):
    seqs = {}
    current_id = None
    current_seq = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    seqs[current_id] = "".join(current_seq)
                header = line[1:].strip()
                seq_id = header.split()[0]
                current_id = seq_id
                current_seq = []
            else:
                current_seq.append(line.strip())

        if current_id is not None:
            seqs[current_id] = "".join(current_seq)

    return seqs


def write_fasta(seqs, out_path, line_width=60):
    with open(out_path, "w") as out:
        for seq_id, seq in seqs.items():
            out.write(f">{seq_id}\n")
            for i in range(0, len(seq), line_width):
                out.write(seq[i:i+line_width] + "\n")


def parse_gtf_attributes(attr_string):
    attrs = {}
    matches = re.finditer(r'(\S+)\s+"([^"]+)"', attr_string)
    for m in matches:
        key = m.group(1)
        value = m.group(2)
        attrs[key] = value
    return attrs


def parse_gtf_transcript_to_gene(gtf_path):
    transcript_to_gene = {}
    gene_to_transcripts = defaultdict(set)

    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = parse_gtf_attributes(fields[8])

            transcript_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id")

            if transcript_id and gene_id:
                transcript_to_gene[transcript_id] = gene_id
                gene_to_transcripts[gene_id].add(transcript_id)

    return transcript_to_gene, gene_to_transcripts


def main():
    parser = argparse.ArgumentParser(
        description="Extract the longest protein isoform per gene from BRAKER outputs."
    )
    parser.add_argument("--gtf", required=True, help="Input BRAKER GTF file")
    parser.add_argument("--aa", required=True, help="Input BRAKER protein FASTA")
    parser.add_argument("--cds", required=True, help="Input BRAKER coding sequence FASTA")
    parser.add_argument("--out-aa", required=True, help="Output longest-isoform protein FASTA")
    parser.add_argument("--out-cds", required=True, help="Output longest-isoform CDS FASTA")
    parser.add_argument("--out-tsv", required=True, help="Output summary TSV")
    args = parser.parse_args()

    print("Parsing GTF...", file=sys.stderr)
    transcript_to_gene, gene_to_transcripts = parse_gtf_transcript_to_gene(args.gtf)

    print("Reading protein FASTA...", file=sys.stderr)
    aa_seqs = parse_fasta(args.aa)

    print("Reading CDS FASTA...", file=sys.stderr)
    cds_seqs = parse_fasta(args.cds)

    missing_in_gtf = []
    protein_lengths = {}
    gene_best = {}

    for transcript_id, seq in aa_seqs.items():
        if transcript_id not in transcript_to_gene:
            missing_in_gtf.append(transcript_id)
            continue

        gene_id = transcript_to_gene[transcript_id]
        prot_len = len(seq.rstrip("*"))  # remove terminal stop if present
        protein_lengths[transcript_id] = prot_len

        if gene_id not in gene_best:
            gene_best[gene_id] = (transcript_id, prot_len)
        else:
            best_tid, best_len = gene_best[gene_id]
            if prot_len > best_len:
                gene_best[gene_id] = (transcript_id, prot_len)
            elif prot_len == best_len:
                # deterministic tie-breaker: keep lexicographically smaller transcript ID
                if transcript_id < best_tid:
                    gene_best[gene_id] = (transcript_id, prot_len)

    selected_transcripts = {tid for tid, _ in gene_best.values()}

    out_aa = {}
    out_cds = {}

    for gene_id, (transcript_id, prot_len) in sorted(gene_best.items()):
        if transcript_id in aa_seqs:
            out_aa[transcript_id] = aa_seqs[transcript_id]
        if transcript_id in cds_seqs:
            out_cds[transcript_id] = cds_seqs[transcript_id]

    print("Writing filtered FASTA files...", file=sys.stderr)
    write_fasta(out_aa, args.out_aa)
    write_fasta(out_cds, args.out_cds)

    print("Writing summary table...", file=sys.stderr)
    with open(args.out_tsv, "w") as out:
        out.write("gene_id\tselected_transcript_id\tprotein_length_aa\tcds_present\tprotein_present\n")
        for gene_id, (transcript_id, prot_len) in sorted(gene_best.items()):
            cds_present = "yes" if transcript_id in cds_seqs else "no"
            protein_present = "yes" if transcript_id in aa_seqs else "no"
            out.write(f"{gene_id}\t{transcript_id}\t{prot_len}\t{cds_present}\t{protein_present}\n")

    total_genes = len(gene_best)
    total_selected = len(selected_transcripts)
    total_aa_written = len(out_aa)
    total_cds_written = len(out_cds)

    print("\nDone.", file=sys.stderr)
    print(f"Genes with selected isoform: {total_genes}", file=sys.stderr)
    print(f"Selected transcripts:        {total_selected}", file=sys.stderr)
    print(f"Protein FASTA written:       {total_aa_written}", file=sys.stderr)
    print(f"CDS FASTA written:           {total_cds_written}", file=sys.stderr)

    if missing_in_gtf:
        print(
            f"Warning: {len(missing_in_gtf)} protein IDs were found in the FASTA but not in the GTF.",
            file=sys.stderr
        )
        print("First few missing IDs:", ", ".join(missing_in_gtf[:10]), file=sys.stderr)


if __name__ == "__main__":
    main()
