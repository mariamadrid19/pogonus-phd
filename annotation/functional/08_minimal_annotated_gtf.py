#!/usr/bin/env python3

import csv
import re
from pathlib import Path

PROJECT = Path("/scratch/leuven/357/vsc35707/annotation/func-annotation")

BRAKER_GTF = PROJECT / "input/braker/braker.gtf"
ANNOT_TSV = PROJECT / "results/06_final_annotation/braker_functional_annotation.compact.tsv"
OUTDIR = PROJECT / "results/08_minimal_annotation"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_GTF = OUTDIR / "braker.minimal.annotated.gtf"
OUT_MAP = OUTDIR / "braker.minimal.annotation_map.tsv"


def parse_gtf_attributes(attr_string, feature=None):
    attr_string = attr_string.strip()
    attrs = {}

    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr_string):
        attrs[m.group(1)] = m.group(2)

    if attrs:
        return attrs

    if feature == "gene" and attr_string:
        attrs["gene_id"] = attr_string
        return attrs

    if feature in ("transcript", "mRNA") and attr_string:
        attrs["transcript_id"] = attr_string
        if ".t" in attr_string:
            attrs["gene_id"] = attr_string.split(".t")[0]
        return attrs

    return attrs


def format_gtf_attributes(attrs):
    preferred = [
        "gene_id",
        "transcript_id",
        "gene_name",
        "product",
        "annotation_confidence",
    ]
    ordered = []
    seen = set()

    for k in preferred:
        if k in attrs and attrs[k] != "":
            ordered.append(k)
            seen.add(k)

    for k in attrs:
        if k not in seen and attrs[k] != "":
            ordered.append(k)

    return " ".join(f'{k} "{attrs[k]}";' for k in ordered)


def clean_product_name(row):
    for key in [
        "best_product_name",
        "swissprot_gene_name",
        "eggnog_preferred_name",
        "interpro_primary_description_short",
    ]:
        val = (row.get(key) or "").strip()
        if val and val.lower() != "uncharacterized protein":
            return val
    return "uncharacterized protein"


annot_by_gene = {}
annot_by_transcript = {}

with open(ANNOT_TSV, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        gid = row["gene_id"]
        tid = row["selected_transcript_id"]
        product = clean_product_name(row)
        conf = (row.get("annotation_confidence") or "").strip()

        annot = {
            "gene_id": gid,
            "transcript_id": tid,
            "gene_name": product,
            "product": product,
            "annotation_confidence": conf,
        }

        annot_by_transcript[tid] = annot
        if gid not in annot_by_gene:
            annot_by_gene[gid] = annot


with open(OUT_MAP, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "gene_id",
        "transcript_id",
        "gene_name",
        "product",
        "annotation_confidence",
    ])
    for gid in sorted(annot_by_gene):
        a = annot_by_gene[gid]
        writer.writerow([
            a["gene_id"],
            a["transcript_id"],
            a["gene_name"],
            a["product"],
            a["annotation_confidence"],
        ])


with open(BRAKER_GTF) as in_fh, open(OUT_GTF, "w") as gtf_out:
    for line in in_fh:
        if not line.strip():
            continue

        if line.startswith("#"):
            gtf_out.write(line)
            continue

        fields = line.rstrip("\n").split("\t")
        if len(fields) != 9:
            continue

        seqid, source, feature, start, end, score, strand, phase, attr_string = fields
        attrs = parse_gtf_attributes(attr_string, feature=feature)

        gene_id = attrs.get("gene_id", "")
        transcript_id = attrs.get("transcript_id", "")

        annot = None
        if transcript_id and transcript_id in annot_by_transcript:
            annot = annot_by_transcript[transcript_id]
        elif gene_id and gene_id in annot_by_gene:
            annot = annot_by_gene[gene_id]

        gtf_attrs = dict(attrs)
        if annot:
            gtf_attrs["gene_name"] = annot["gene_name"]
            gtf_attrs["product"] = annot["product"]
            gtf_attrs["annotation_confidence"] = annot["annotation_confidence"]

        gtf_out.write("\t".join([
            seqid, source, feature, start, end, score, strand, phase,
            format_gtf_attributes(gtf_attrs)
        ]) + "\n")

print(f"Wrote minimal annotated GTF: {OUT_GTF}")
print(f"Wrote minimal annotation map: {OUT_MAP}")
print(f"Annotated transcripts: {len(annot_by_transcript)}")
print(f"Annotated genes: {len(annot_by_gene)}")
