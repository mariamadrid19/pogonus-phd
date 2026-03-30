#!/usr/bin/env python3

import csv
import re
from pathlib import Path

PROJECT = Path("/scratch/leuven/357/vsc35707/annotation/func-annotation")

BRAKER_GTF = PROJECT / "input/braker/braker.gtf"
ANNOT_TSV = PROJECT / "results/06_final_annotation/braker_functional_annotation.compact.tsv"
OUTDIR = PROJECT / "results/07_annotated_annotation"
OUTDIR.mkdir(parents=True, exist_ok=True)

OUT_GTF = OUTDIR / "braker.annotated.gtf"
OUT_GFF3 = OUTDIR / "braker.annotated.gff3"
OUT_MAP = OUTDIR / "braker.annotation_map.tsv"


def parse_gtf_attributes(attr_string, feature=None):
    attr_string = attr_string.strip()
    attrs = {}

    # Standard GTF attributes
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr_string):
        attrs[m.group(1)] = m.group(2)

    if attrs:
        return attrs

    # BRAKER shorthand:
    # gene line:       g1
    # transcript line: g1.t1
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
        "note",
        "annotation_confidence",
    ]
    ordered = []
    seen = set()

    for k in preferred:
        if k in attrs:
            ordered.append(k)
            seen.add(k)

    for k in attrs:
        if k not in seen:
            ordered.append(k)

    return " ".join(f'{k} "{attrs[k]}";' for k in ordered if attrs[k] != "")


def gff3_escape(text):
    if text is None:
        return ""
    text = str(text)
    text = text.replace("%", "%25")
    text = text.replace(";", "%3B")
    text = text.replace("=", "%3D")
    text = text.replace("&", "%26")
    text = text.replace(",", "%2C")
    text = text.replace("\t", " ")
    return text


def format_gff3_attributes(attrs):
    return ";".join(
        f"{k}={gff3_escape(v)}" for k, v in attrs.items() if v not in (None, "")
    )


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


def build_note(row):
    parts = []

    desc = (row.get("eggnog_description_short") or "").strip()
    if desc:
        parts.append(f"eggNOG: {desc}")

    ipr = (row.get("interpro_primary_id") or "").strip()
    ipr_desc = (row.get("interpro_primary_description_short") or "").strip()
    if ipr and ipr_desc:
        parts.append(f"InterPro: {ipr} {ipr_desc}")
    elif ipr:
        parts.append(f"InterPro: {ipr}")

    swiss = (row.get("swissprot_stitle_short") or "").strip()
    if swiss:
        parts.append(f"Swiss-Prot: {swiss}")

    source = (row.get("best_product_source") or "").strip()
    if source:
        parts.append(f"best_source={source}")

    return " | ".join(parts)


# Read compact annotation and build gene-level and transcript-level maps
annot_by_gene = {}
annot_by_transcript = {}

with open(ANNOT_TSV, newline="") as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    for row in reader:
        gid = row["gene_id"]
        tid = row["selected_transcript_id"]

        annot = {
            "gene_id": gid,
            "transcript_id": tid,
            "gene_name": clean_product_name(row),
            "product": clean_product_name(row),
            "note": build_note(row),
            "annotation_confidence": (row.get("annotation_confidence") or "").strip(),
        }

        annot_by_transcript[tid] = annot

        # Keep one best annotation per gene
        if gid not in annot_by_gene:
            annot_by_gene[gid] = annot


# Write annotation map
with open(OUT_MAP, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow([
        "gene_id",
        "transcript_id",
        "gene_name",
        "product",
        "note",
        "annotation_confidence"
    ])
    for gid in sorted(annot_by_gene):
        a = annot_by_gene[gid]
        writer.writerow([
            a["gene_id"],
            a["transcript_id"],
            a["gene_name"],
            a["product"],
            a["note"],
            a["annotation_confidence"],
        ])


with open(BRAKER_GTF) as in_fh, open(OUT_GTF, "w") as gtf_out, open(OUT_GFF3, "w") as gff3_out:
    gff3_out.write("##gff-version 3\n")

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

        # Prefer transcript-level annotation if exact transcript exists,
        # otherwise fall back to gene-level annotation
        annot = None
        if transcript_id and transcript_id in annot_by_transcript:
            annot = annot_by_transcript[transcript_id]
        elif gene_id and gene_id in annot_by_gene:
            annot = annot_by_gene[gene_id]

        # ---------- GTF ----------
        gtf_attrs = dict(attrs)

        if annot:
            gtf_attrs["gene_name"] = annot["gene_name"]
            gtf_attrs["product"] = annot["product"]
            gtf_attrs["note"] = annot["note"]
            gtf_attrs["annotation_confidence"] = annot["annotation_confidence"]

        gtf_out.write("\t".join([
            seqid, source, feature, start, end, score, strand, phase,
            format_gtf_attributes(gtf_attrs)
        ]) + "\n")

        # ---------- GFF3 ----------
        gff3_feature = feature
        gff3_attrs = {}

        if feature == "gene":
            gene_feature_id = gene_id if gene_id else f"gene:{seqid}:{start}-{end}"
            gff3_attrs["ID"] = gene_feature_id
            if annot:
                gff3_attrs["Name"] = annot["gene_name"]
                gff3_attrs["gene_name"] = annot["gene_name"]
                gff3_attrs["product"] = annot["product"]
                gff3_attrs["note"] = annot["note"]
                gff3_attrs["annotation_confidence"] = annot["annotation_confidence"]

        elif feature in ("transcript", "mRNA"):
            gff3_feature = "mRNA"
            tx_feature_id = transcript_id if transcript_id else f"transcript:{seqid}:{start}-{end}"
            gff3_attrs["ID"] = tx_feature_id
            if gene_id:
                gff3_attrs["Parent"] = gene_id
            if annot:
                gff3_attrs["Name"] = annot["gene_name"]
                gff3_attrs["gene_name"] = annot["gene_name"]
                gff3_attrs["product"] = annot["product"]
                gff3_attrs["note"] = annot["note"]
                gff3_attrs["annotation_confidence"] = annot["annotation_confidence"]

        else:
            parent = transcript_id if transcript_id else gene_id
            if parent:
                gff3_attrs["ID"] = f"{feature}:{parent}:{start}-{end}"
                gff3_attrs["Parent"] = parent
            if annot and feature in ("CDS", "exon", "start_codon", "stop_codon"):
                gff3_attrs["gene_name"] = annot["gene_name"]
                gff3_attrs["product"] = annot["product"]
                gff3_attrs["annotation_confidence"] = annot["annotation_confidence"]

        gff3_out.write("\t".join([
            seqid, source, gff3_feature, start, end, score, strand, phase,
            format_gff3_attributes(gff3_attrs)
        ]) + "\n")

print(f"Wrote annotated GTF:  {OUT_GTF}")
print(f"Wrote annotated GFF3: {OUT_GFF3}")
print(f"Wrote map table:      {OUT_MAP}")
print(f"Annotated transcripts: {len(annot_by_transcript)}")
print(f"Annotated genes:       {len(annot_by_gene)}")
