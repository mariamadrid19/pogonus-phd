#!/usr/bin/env python3

import csv
from pathlib import Path
from collections import defaultdict

PROJECT = Path("/scratch/leuven/357/vsc35707/annotation/func-annotation")

LONGEST_TSV = PROJECT / "results/01_longest_isoforms/braker.longest_isoforms.tsv"
EGGNOG_TSV = PROJECT / "results/03_eggnog/braker_longest_isoforms.emapper.annotations"
INTERPRO_TSV = PROJECT / "results/04_interproscan/braker_longest_isoforms.tsv"
DIAMOND_TSV = PROJECT / "results/05_diamond/braker_longest_isoforms_vs_swissprot.tsv"

OUTDIR = PROJECT / "results/06_final_annotation"
OUTDIR.mkdir(parents=True, exist_ok=True)

FULL_OUT = OUTDIR / "braker_functional_annotation.full.tsv"
COMPACT_OUT = OUTDIR / "braker_functional_annotation.compact.tsv"
STATS_OUT = OUTDIR / "braker_functional_annotation.summary.txt"


def clean_value(v: str) -> str:
    if v is None:
        return ""
    v = v.strip()
    if v == "-" or v == "":
        return ""
    return v


def uniq_preserve(values):
    seen = set()
    out = []
    for v in values:
        v = clean_value(v)
        if not v:
            continue
        if v not in seen:
            seen.add(v)
            out.append(v)
    return out


def join_unique(values, sep=";"):
    vals = uniq_preserve(values)
    return sep.join(vals) if vals else ""


def truncate_text(text: str, max_len: int = 300) -> str:
    text = clean_value(text)
    if len(text) <= max_len:
        return text
    return text[: max_len - 3] + "..."


def parse_longest_isoforms(path):
    data = {}
    with open(path, newline="") as fh:
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
    header = None

    with open(path, newline="") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").rstrip("\n").split("\t")
                break

        if header is None:
            raise RuntimeError("Could not find eggNOG header line.")

        reader = csv.DictReader(fh, fieldnames=header, delimiter="\t")

        for row in reader:
            tid = row.get("query", "").strip()
            if not tid:
                continue

            data[tid] = {
                "eggnog_seed_ortholog": clean_value(row.get("seed_ortholog", "")),
                "eggnog_evalue": clean_value(row.get("evalue", "")),
                "eggnog_score": clean_value(row.get("score", "")),
                "eggnog_ogs": clean_value(row.get("eggNOG_OGs", "")),
                "eggnog_max_annot_lvl": clean_value(row.get("max_annot_lvl", "")),
                "eggnog_cog_category": clean_value(row.get("COG_category", "")),
                "eggnog_preferred_name": clean_value(row.get("Preferred_name", "")),
                "eggnog_description": clean_value(row.get("Description", "")),
                "eggnog_go": clean_value(row.get("GOs", "")),
                "eggnog_ec": clean_value(row.get("EC", "")),
                "eggnog_kegg_ko": clean_value(row.get("KEGG_ko", "")),
                "eggnog_kegg_pathway": clean_value(row.get("KEGG_Pathway", "")),
                "eggnog_kegg_module": clean_value(row.get("KEGG_Module", "")),
                "eggnog_kegg_reaction": clean_value(row.get("KEGG_Reaction", "")),
                "eggnog_kegg_rclass": clean_value(row.get("KEGG_rclass", "")),
                "eggnog_brite": clean_value(row.get("BRITE", "")),
                "eggnog_kegg_tc": clean_value(row.get("KEGG_TC", "")),
                "eggnog_cazy": clean_value(row.get("CAZy", "")),
                "eggnog_bigg_reaction": clean_value(row.get("BiGG_Reaction", "")),
                "eggnog_pfam": clean_value(row.get("PFAMs", "")),
            }

    return data


def parse_interpro(path):
    analyses = defaultdict(list)
    accessions = defaultdict(list)
    descriptions = defaultdict(list)
    go_terms = defaultdict(list)
    pathways = defaultdict(list)

    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue

            tid = row[0]
            analysis = row[3] if len(row) > 3 else ""
            ipr_id = row[11] if len(row) > 11 else ""
            ipr_desc = row[12] if len(row) > 12 else ""
            go_field = row[13] if len(row) > 13 else ""
            pathway_field = row[14] if len(row) > 14 else ""

            if clean_value(analysis):
                analyses[tid].append(analysis)
            if clean_value(ipr_id):
                accessions[tid].append(ipr_id)
            if clean_value(ipr_desc):
                descriptions[tid].append(ipr_desc)

            if clean_value(go_field):
                go_terms[tid].extend([x for x in go_field.split("|") if clean_value(x)])
            if clean_value(pathway_field):
                pathways[tid].extend([x for x in pathway_field.split("|") if clean_value(x)])

    data = {}
    tids = set(analyses) | set(accessions) | set(descriptions) | set(go_terms) | set(pathways)

    for tid in tids:
        uniq_desc = uniq_preserve(descriptions[tid])
        uniq_ipr = uniq_preserve(accessions[tid])
        uniq_go = uniq_preserve(go_terms[tid])
        uniq_path = uniq_preserve(pathways[tid])
        uniq_analyses = uniq_preserve(analyses[tid])

        data[tid] = {
            "interpro_analyses": ";".join(uniq_analyses),
            "interpro_ids": ";".join(uniq_ipr),
            "interpro_descriptions": ";".join(uniq_desc),
            "interpro_go": ";".join(uniq_go),
            "interpro_pathways": ";".join(uniq_path),

            # Compact-friendly summaries
            "interpro_primary_id": uniq_ipr[0] if uniq_ipr else "",
            "interpro_primary_description": uniq_desc[0] if uniq_desc else "",
            "interpro_id_count": str(len(uniq_ipr)) if uniq_ipr else "0",
            "interpro_go_count": str(len(uniq_go)) if uniq_go else "0",
            "interpro_pathway_count": str(len(uniq_path)) if uniq_path else "0",
        }

    return data


def parse_diamond(path):
    data = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row:
                continue

            tid = row[0]
            if tid in data:
                continue  # keep best hit only

            stitle = row[12] if len(row) > 12 else ""
            subject = row[1]

            swissprot_gene_name = ""
            if " GN=" in stitle:
                swissprot_gene_name = stitle.split(" GN=", 1)[1].split()[0]

            data[tid] = {
                "swissprot_subject": clean_value(subject),
                "swissprot_pident": clean_value(row[2]),
                "swissprot_align_length": clean_value(row[3]),
                "swissprot_mismatch": clean_value(row[4]),
                "swissprot_gapopen": clean_value(row[5]),
                "swissprot_qstart": clean_value(row[6]),
                "swissprot_qend": clean_value(row[7]),
                "swissprot_sstart": clean_value(row[8]),
                "swissprot_send": clean_value(row[9]),
                "swissprot_evalue": clean_value(row[10]),
                "swissprot_bitscore": clean_value(row[11]),
                "swissprot_stitle": clean_value(stitle),
                "swissprot_gene_name": swissprot_gene_name,
            }

    return data


def assign_confidence(row):
    has_swiss = bool(clean_value(row.get("swissprot_subject", "")))
    has_ipr = bool(clean_value(row.get("interpro_ids", "")))
    has_eggnog = bool(
        clean_value(row.get("eggnog_preferred_name", "")) or
        clean_value(row.get("eggnog_description", ""))
    )

    if has_swiss and has_ipr:
        return "high"
    if has_swiss or (has_ipr and has_eggnog):
        return "medium"
    if has_ipr or has_eggnog:
        return "low"
    return "uncharacterized"


def best_product_name(row):
    swiss_gene = clean_value(row.get("swissprot_gene_name", ""))
    eggnog_name = clean_value(row.get("eggnog_preferred_name", ""))
    interpro_desc = clean_value(row.get("interpro_primary_description", ""))

    if swiss_gene:
        return swiss_gene
    if eggnog_name:
        return eggnog_name
    if interpro_desc:
        return interpro_desc
    return "uncharacterized protein"


def best_product_source(row):
    if clean_value(row.get("swissprot_gene_name", "")):
        return "Swiss-Prot"
    if clean_value(row.get("eggnog_preferred_name", "")):
        return "eggNOG"
    if clean_value(row.get("interpro_primary_description", "")):
        return "InterPro"
    return "none"


longest = parse_longest_isoforms(LONGEST_TSV)
eggnog = parse_eggnog(EGGNOG_TSV)
interpro = parse_interpro(INTERPRO_TSV)
diamond = parse_diamond(DIAMOND_TSV)

full_fields = [
    "gene_id",
    "selected_transcript_id",
    "protein_length_aa",
    "cds_present",
    "protein_present",
    "best_product_name",
    "best_product_source",
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
    "interpro_primary_id",
    "interpro_primary_description",
    "interpro_id_count",
    "interpro_go_count",
    "interpro_pathway_count",
    "swissprot_subject",
    "swissprot_gene_name",
    "swissprot_pident",
    "swissprot_align_length",
    "swissprot_evalue",
    "swissprot_bitscore",
    "swissprot_stitle",
    "annotation_confidence",
]

compact_fields = [
    "gene_id",
    "selected_transcript_id",
    "protein_length_aa",
    "best_product_name",
    "best_product_source",
    "eggnog_preferred_name",
    "eggnog_description_short",
    "interpro_primary_id",
    "interpro_primary_description_short",
    "interpro_id_count",
    "swissprot_subject",
    "swissprot_gene_name",
    "swissprot_pident",
    "swissprot_evalue",
    "swissprot_bitscore",
    "swissprot_stitle_short",
    "annotation_confidence",
]

rows_full = []
rows_compact = []

for tid in sorted(longest):
    row = {}
    row.update(longest.get(tid, {}))
    row.update(eggnog.get(tid, {}))
    row.update(interpro.get(tid, {}))
    row.update(diamond.get(tid, {}))

    row["annotation_confidence"] = assign_confidence(row)
    row["best_product_name"] = best_product_name(row)
    row["best_product_source"] = best_product_source(row)

    row["eggnog_description_short"] = truncate_text(row.get("eggnog_description", ""), 300)
    row["interpro_primary_description_short"] = truncate_text(row.get("interpro_primary_description", ""), 200)
    row["swissprot_stitle_short"] = truncate_text(row.get("swissprot_stitle", ""), 200)

    rows_full.append({k: row.get(k, "") for k in full_fields})
    rows_compact.append({k: row.get(k, "") for k in compact_fields})

with open(FULL_OUT, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=full_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows_full)

with open(COMPACT_OUT, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=compact_fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(rows_compact)

n_total = len(longest)
n_eggnog = sum(1 for t in longest if t in eggnog)
n_interpro = sum(1 for t in longest if t in interpro)
n_swiss = sum(1 for t in longest if t in diamond)

conf_counts = defaultdict(int)
for row in rows_full:
    conf_counts[row["annotation_confidence"]] += 1

with open(STATS_OUT, "w") as fh:
    fh.write(f"Representative proteins\t{n_total}\n")
    fh.write(f"With eggNOG annotation\t{n_eggnog}\n")
    fh.write(f"With InterPro annotation\t{n_interpro}\n")
    fh.write(f"With Swiss-Prot hit\t{n_swiss}\n")
    fh.write(f"High confidence\t{conf_counts['high']}\n")
    fh.write(f"Medium confidence\t{conf_counts['medium']}\n")
    fh.write(f"Low confidence\t{conf_counts['low']}\n")
    fh.write(f"Uncharacterized\t{conf_counts['uncharacterized']}\n")

print(f"Wrote full table:    {FULL_OUT}")
print(f"Wrote compact table: {COMPACT_OUT}")
print(f"Wrote summary:       {STATS_OUT}")
print(f"Representative proteins: {n_total}")
print(f"With eggNOG annotation:  {n_eggnog}")
print(f"With InterPro annotation:{n_interpro}")
print(f"With Swiss-Prot hit:     {n_swiss}")
print(
    "Confidence counts:",
    f"high={conf_counts['high']},",
    f"medium={conf_counts['medium']},",
    f"low={conf_counts['low']},",
    f"uncharacterized={conf_counts['uncharacterized']}"
)
