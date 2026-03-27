# Functional annotation pipeline for BRAKER gene predictions

This repository contains a step-by-step workflow to functionally annotate a genome starting from **BRAKER** output files.

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
   Build a final annotation table that can be joined back to the genome annotation.

---

## Input files

The pipeline starts from BRAKER output:

- `braker.aa`  
  Predicted protein sequences. This is the main input for functional annotation.

- `braker.codingseq`  
  Predicted CDS nucleotide sequences. Useful for keeping a CDS file matched to the filtered protein set.

- `braker.gtf`  
  Structural annotation with gene, transcript, exon, and CDS coordinates.
