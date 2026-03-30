#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

PROJECT = Path("/scratch/leuven/357/vsc35707/annotation/func-annotation")

LONGEST_TSV = PROJECT / "results/01_longest_isoforms/braker.longest_isoforms.tsv"
EGGNOG_TSV = PROJECT / "results/03_eggnog/braker_longest_isoforms.emapper.annotations"
INTERPRO_TSV = PROJECT / "results/04_interproscan/braker_longest_isoforms.tsv"
DIAMOND_TSV = PROJECT / "results/05_diamond/braker_longest_isoforms_vs_swissprot.tsv"
OUTDIR = PROJECT / "results/06_final_annotation"
OUTFILE = OUTDIR / "braker_functional_annotation.tsv"

OUTDIR.mkdir(parents=True, exist_ok=True)


def uniq_join(values):
    cleaned = []
    seen = set()
    for v in values:
        if v is None:
            continue
        v = v.strip()
        if not v or v == "-" or v in seen:
            continue
        seen.add(v)
        cleaned.append(v)
    return ";".join(cleaned) if cleaned else ""


def parse_longest_isoforms(path):
    data = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            tid = row["selected_transcript_id"]
            data[tid] = {
                "gene_id": row["gene_id"],
                "selected_transcript_id": tid,
                "protein_length_aa": row["protein_length_aa"],
                "cds_present": row["cds_present"],
                "protein_present": row["protein_present"],
            }
    return data


def parse_eggnog(path):
    data = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break
        reader = csv.DictReader(fh, fieldnames=header, delimiter="\t")
        for row in reader:
            tid = row["query"]
            data[tid] = {
                "eggnog_seed_ortholog": row.get("seed_ortholog", ""),
                "eggnog_evalue": row.get("evalue", ""),
                "eggnog_score": row.get("score", ""),
                "eggnog_ogs": row.get("eggNOG_OGs", ""),
                "eggnog_max_annot_lvl": row.get("max_annot_lvl", ""),
                "eggnog_cog_category": row.get("COG_category", ""),
                "eggnog_description": row.get("Description", ""),
                "eggnog_preferred_name": row.get("Preferred_name", ""),
                "eggnog_go": row.get("GOs", "") if "GOs" in row else row.get("GOsEC", ""),
                "eggnog_ec": row.get("EC", "") if "EC" in row else row.get("GOsEC", ""),
                "eggnog_kegg_ko": row.get("KEGG_ko", ""),
                "eggnog_kegg_pathway": row.get("KEGG_Pathway", ""),
                "eggnog_kegg_module": row.get("KEGG_Module", ""),
                "eggnog_kegg_reaction": row.get("KEGG_Reaction", ""),
                "eggnog_kegg_rclass": row.get("KEGG_rclass", ""),
                "eggnog_brite": row.get("BRITE", ""),
                "eggnog_kegg_tc": row.get("KEGG_TC", ""),
                "eggnog_cazy": row.get("CAZy", ""),
                "eggnog_bigg_reaction": row.get("BiGG_Reaction", ""),
                "eggnog_pfam": row.get("PFAMs", ""),
            }
    return data


def parse_interpro(path):
    ids = defaultdict(list)
    descs = defaultdict(list)
    gos = defaultdict(list)
    paths = defaultdict(list)
    analyses = defaultdict(list)

    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            tid = row[0]
            analysis = row[3] if len(row) > 3 else ""
            ipr_id = row[11] if len(row) > 11 else ""
            ipr_desc = row[12] if len(row) > 12 else ""
            go_terms = row[13] if len(row) > 13 else ""
            pathway = row[14] if len(row) > 14 else ""

            analyses[tid].append(analysis)
            if ipr_id and ipr_id != "-":
                ids[tid].append(ipr_id)
            if ipr_desc and ipr_desc != "-":
                descs[tid].append(ipr_desc)
            if go_terms and go_terms != "-":
                gos[tid].extend(go_terms.split("|"))
            if pathway and pathway != "-":
                paths[tid].extend(pathway.split("|"))

    data = {}
    all_tids = set(analyses) | set(ids) | set(descs) | set(gos) | set(paths)
    for tid in all_tids:
        data[tid] = {
            "interpro_analyses": uniq_join(analyses[tid]),
            "interpro_ids": uniq_join(ids[tid]),
            "interpro_descriptions": uniq_join(descs[tid]),
            "interpro_go": uniq_join(gos[tid]),
            "interpro_pathways": uniq_join(paths[tid]),
        }
    return data


def parse_diamond(path):
    data = {}
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue
            tid = row[0]
            if tid in data:
                continue
            data[tid] = {
                "swissprot_subject": row[1],
                "swissprot_pident": row[2],
                "swissprot_align_length": row[3],
                "swissprot_mismatch": row[4],
                "swissprot_gapopen": row[5],
                "swissprot_qstart": row[6],
                "swissprot_qend": row[7],
                "swissprot_sstart": row[8],
                "swissprot_send": row[9],
                "swissprot_evalue": row[10],
                "swissprot_bitscore": row[11],
                "swissprot_stitle": row[12] if len(row) > 12 else "",
            }
    return data


def assign_confidence(row):
    has_swiss = bool(row.get("swissprot_subject"))
    has_ipr = bool(row.get("interpro_ids"))
    has_eggnog = bool(row.get("eggnog_preferred_name") or row.get("eggnog_description"))

    if has_swiss and has_ipr:
        return "high"
    if has_swiss or (has_ipr and has_eggnog):
        return "medium"
    if has_ipr or has_eggnog:
        return "low"
    return "uncharacterized"


longest = parse_longest_isoforms(LONGEST_TSV)
eggnog = parse_eggnog(EGGNOG_TSV)
interpro = parse_interpro(INTERPRO_TSV)
diamond = parse_diamond(DIAMOND_TSV)

fieldnames = [
    "gene_id",
    "selected_transcript_id",
    "protein_length_aa",
    "cds_present",
    "protein_present",
    "eggnog_seed_ortholog",
    "eggnog_evalue",
    "eggnog_score",
    "eggnog_ogs",
    "eggnog_max_annot_lvl",
    "eggnog_cog_category",
    "eggnog_preferred_name",
    "eggnog_description",
    "eggnog_go",
    "eggnog_ec",
    "eggnog_kegg_ko",
    "eggnog_kegg_pathway",
    "eggnog_kegg_module",
    "eggnog_kegg_reaction",
    "eggnog_kegg_rclass",
    "eggnog_brite",
    "eggnog_kegg_tc",
    "eggnog_cazy",
    "eggnog_bigg_reaction",
    "eggnog_pfam",
    "interpro_analyses",
    "interpro_ids",
    "interpro_descriptions",
    "interpro_go",
    "interpro_pathways",
    "swissprot_subject",
    "swissprot_pident",
    "swissprot_align_length",
    "swissprot_evalue",
    "swissprot_bitscore",
    "swissprot_stitle",
    "annotation_confidence",
]

with open(OUTFILE, "w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()

    for tid in sorted(longest):
        row = {}
        row.update(longest.get(tid, {}))
        row.update(eggnog.get(tid, {}))
        row.update(interpro.get(tid, {}))
        row.update(diamond.get(tid, {}))
        row["annotation_confidence"] = assign_confidence(row)
        writer.writerow({k: row.get(k, "") for k in fieldnames})

print(f"Wrote: {OUTFILE}")
print(f"Representative proteins: {len(longest)}")
print(f"With eggNOG annotation: {sum(1 for t in longest if t in eggnog)}")
print(f"With InterPro annotation: {sum(1 for t in longest if t in interpro)}")
print(f"With Swiss-Prot hit: {sum(1 for t in longest if t in diamond)}")
