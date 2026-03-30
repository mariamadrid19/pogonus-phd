# Functional annotation pipeline for BRAKER3 gene predictions

This repository contains a step-by-step workflow to functionally annotate a genome starting from **BRAKER3** output files.

The goal is to go from raw structural predictions such as:

- `braker.aa`
- `braker.codingseq`
- `braker.gtf`

to a final annotation table describing what the predicted genes are likely to encode.

---

## Overview

BRAKER predicts **gene models**. That tells us **where genes are** in the genome and what their predicted transcripts and proteins look like.

Functional annotation is the next layer: it helps answer questions like:

- What proteins do these genes encode?
- Which genes have known domains?
- Which genes have orthologs in other organisms?
- Which genes can be linked to GO terms, pathways, or enzyme functions?

This pipeline is built around the idea that the most reliable annotation comes from combining multiple sources of evidence:

1. **Representative protein selection**  
   Keep one protein isoform per gene to avoid inflating downstream results.

2. **Gene set quality assessment**  
   Use BUSCO to check whether the predicted proteome is reasonably complete.

3. **Orthology-based annotation**  
   Use eggNOG-mapper to assign broad functional descriptions, GO terms, and pathway terms.

4. **Domain-based annotation**  
   Use InterProScan to identify conserved domains and protein families.

5. **Sequence similarity against curated proteins**  
   Use DIAMOND against Swiss-Prot for high-confidence naming and manual curation.

6. **Merge all evidence**  
   Build final annotation tables that can be joined back to the genome annotation or used directly in downstream analyses.

7. **Decorate the BRAKER annotation**  
   Add functional labels back onto the original annotation to generate richer GTF/GFF3 files.

8. **Create a minimal downstream-friendly GTF**  
   Produce a lighter GTF with only the most useful attributes for RNA-seq and related analyses.

---

## Input files

The pipeline starts from BRAKER output:

- `braker.aa`  
  Predicted protein sequences. This is the main input for functional annotation.

- `braker.codingseq`  
  Predicted CDS nucleotide sequences. Useful for keeping a CDS file matched to the filtered protein set.

- `braker.gtf`  
  Structural annotation with gene, transcript, exon, and CDS coordinates.

---

## Scripts in this repository

### `01a_longest_isoform.py`

This script selects one representative transcript per gene, using the **longest protein isoform** as the default representative. It parses the BRAKER GTF to connect transcripts to genes, measures protein length from the protein FASTA, and keeps a single protein per gene.

This step reduces redundancy before downstream annotation. Many annotation tools work better on one representative isoform per gene, and this also makes later summary tables easier to interpret.

The script writes:
- a filtered protein FASTA
- a filtered CDS FASTA
- a summary table showing which transcript was selected for each gene

---

### `01b_longest_isoform.sh`

This is the SLURM wrapper used to run the longest-isoform extraction step on the cluster.

It exists to make the first step reproducible and easy to rerun in batch mode, while writing a log of the run.

---

### `02_busco.sh`

This script runs **BUSCO** on the representative protein set.

BUSCO assesses how complete the predicted proteome is by searching for highly conserved single-copy orthologs expected in a chosen lineage. This is used as a quality check before moving on to functional annotation.

A strong BUSCO result gives confidence that the gene set is suitable for downstream annotation.

---

### `03a_download_eggnog.sh`

This script downloads the required **eggNOG-mapper databases**.

eggNOG-mapper needs its local annotation database, taxonomy database, and DIAMOND database before it can run. This script retrieves those files and organizes them locally for the annotation step.

This database setup is kept separate from the actual annotation run so it only needs to be done once.

---

### `03b_run_eggnog.sh`

This script runs **eggNOG-mapper** on the representative protein FASTA.

eggNOG-mapper performs **orthology-based functional annotation**, providing putative gene names, descriptions, GO terms, KEGG terms, and other broad functional categories. This is often the fastest way to get a first-pass annotation of a newly predicted genome.

This step produces the main orthology-based annotation table used later in the merge step.

---

### `04_interproscan.sh`

This script runs **InterProScan** on the representative proteins.

InterProScan identifies conserved domains, motifs, families, and functional signatures across many member databases. It provides an independent source of evidence that complements eggNOG.

Because InterProScan does not accept `*` characters in protein sequences, the workflow first generates a cleaned protein FASTA with stop symbols removed and then uses that cleaned file as input.

This step is especially useful for:
- confirming functional assignments
- identifying domain architecture
- assigning conservative family-level labels when no precise name is available

---

### `05_diamond.sh`

This script runs **DIAMOND blastp** against **UniProtKB/Swiss-Prot**.

Swiss-Prot is a curated protein database, so this step is mainly used for more reliable naming and manual validation. It helps support or refine the annotation suggested by eggNOG and InterPro.

This step produces a tabular file of best similarity hits that is later merged into the final annotation table.

---

### `06a_merge_annotations.py`

This script merges all annotation evidence into unified output tables.

It combines:
- the longest-isoform summary
- eggNOG annotations
- InterProScan results
- DIAMOND Swiss-Prot hits

The script collapses the InterProScan output to one row per representative transcript and then builds merged annotation tables.

It produces:
- a **full annotation table** with detailed fields from all tools
- a **compact annotation table** with the most useful summary columns
- a **summary file** with counts of annotated proteins and confidence categories

This is the main integration step of the functional annotation workflow.

---

### `06b_merge.sh`

This is the SLURM wrapper used to run the merge script on the cluster.

It makes the merge step reproducible and keeps a separate log of the final table generation.

---

### `07_make_annotated_gtf.py`

This script takes the merged functional annotation and adds it back onto the original BRAKER annotation.

It produces:
- an **annotated GTF**
- an **annotated GFF3**
- an **annotation map table**

The script is designed to handle the way BRAKER writes `gene` and `transcript` features in the GTF, including shorthand identifiers such as `g1` and `g1.t1`. It assigns one best functional label per gene and propagates that information across all related features.

This step is useful when a coordinate-based annotation file with embedded function labels is needed for browsing, filtering, or interpretation.

---

### `08_minimal_annotated_gtf.py`

This script creates a lighter, downstream-friendly GTF for RNA-seq and related analyses.

Instead of carrying long note fields, it adds only a minimal set of useful attributes:

- `gene_id`
- `transcript_id`
- `gene_name`
- `product`
- `annotation_confidence`

This produces a cleaner annotation file that is easier to use with downstream tools while still preserving the most informative functional labels.

It also writes a small annotation map table that can be used later to annotate differential expression results or other gene-level outputs.

---

## Main outputs

By the end of the pipeline, the main outputs are:

- a representative protein set
- a BUSCO quality assessment
- eggNOG functional annotations
- InterProScan domain annotations
- Swiss-Prot similarity hits
- merged annotation tables
- an annotated GTF
- an annotated GFF3
- a minimal annotated GTF for downstream analyses

---

## Why this workflow is useful

This workflow combines multiple complementary lines of evidence rather than relying on a single annotation tool. That makes the final annotation more robust and easier to interpret.

In particular:

- **eggNOG-mapper** provides broad orthology-based function transfer
- **InterProScan** provides domain-based support
- **Swiss-Prot** provides curated similarity-based naming
- the **merge step** integrates them into a single usable resource
- the **annotated GTF/GFF3** and **minimal GTF** make the results easier to use in downstream genomics and transcriptomics analyses

---

## Recommended usage

A practical way to use the final outputs is:

- use the **full annotation table** for in-depth exploration
- use the **compact annotation table** for quick inspection and summaries
- use the **annotated GTF/GFF3** when a coordinate-based annotation with function labels is needed
- use the **minimal annotated GTF** for downstream RNA-seq analyses
- use the **minimal annotation map** to join functional labels onto DEG or gene-level result tables

---

## Summary

This repository implements a complete functional annotation workflow starting from BRAKER structural predictions and ending with integrated, interpretable gene annotations.

The pipeline is designed to be:

- modular
- reproducible
- cluster-friendly
- practical for downstream analyses

It can be used both to generate a rich functional annotation resource and to produce lightweight annotation files suitable for routine RNA-seq and gene-level analyses.
