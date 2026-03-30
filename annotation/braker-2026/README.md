# BRAKER3 genome annotation pipeline

This repository contains a step-by-step workflow to generate a **structural genome annotation** using repeat masking, RNA-seq evidence, Iso-Seq evidence, and **BRAKER**.

The goal is to go from a genome assembly plus transcriptomic evidence to predicted gene models such as:

- `braker.gtf`
- `braker.aa`
- `braker.codingseq`

These outputs describe where genes are in the genome, how they are structured, and what proteins they are predicted to encode.

---

## Overview

BRAKER is a genome annotation pipeline that predicts **gene models** by combining ab initio prediction with external evidence, especially RNA alignments. In practice, a good BRAKER run depends heavily on the quality of the genome assembly, the repeat masking, and the transcript evidence provided to it.

This pipeline is built around the idea that structural annotation works best when evidence is prepared carefully in separate steps:

1. **Soft-mask the genome**  
   Mask repetitive regions so gene predictors are less likely to mistake repeats for real genes.

2. **Prepare RNA-seq splice-aware alignments**  
   Map short-read RNA-seq data to the genome to provide exon and splice-junction evidence.

3. **Map Iso-Seq long reads**  
   Align long transcript reads to improve transcript structure evidence and support exon connectivity.

4. **Quality-control the alignments**  
   Check that BAM files are valid, sorted, indexed, and suitable for downstream use.

5. **Merge RNA evidence**  
   Combine alignments into a consolidated evidence set for annotation.

6. **Run BRAKER**  
   Use the masked genome and transcript alignments to predict gene structures and generate the final annotation files.

---

## Input files

The pipeline starts from a genome assembly and transcriptomic evidence.

Typical key inputs are:

- a genome FASTA file  
- short-read RNA-seq data  
- Iso-Seq long-read data  

The final goal is to produce BRAKER annotation outputs including:

- `braker.gtf`
- `braker.aa`
- `braker.codingseq`

---

## Pipeline steps

### `01_softmask.sh`

This script soft-masks the genome assembly. Soft-masking converts repetitive sequence into lowercase rather than fully removing it, which allows aligners and annotation tools to still access the sequence while signaling that those regions are repetitive.

This step is done because repeats can confuse gene prediction. Without masking, transposable elements and other repetitive regions may be misannotated as genes. Soft-masking is therefore an important preprocessing step that improves the specificity of structural annotation while preserving the underlying sequence information.

### `02a_star_index.sh`

This script builds a STAR genome index for the assembly. STAR requires a precomputed index before RNA-seq reads can be aligned.

This step is done because RNA-seq evidence is one of the main sources BRAKER uses to improve gene prediction. Creating the genome index is the preparation step that makes splice-aware mapping possible in the next stage. It converts the reference genome into a format optimized for rapid RNA-seq alignment.

### `02b_map_star.sh`

This script maps short-read RNA-seq data to the genome using STAR. Because STAR is splice-aware, it can align reads across exon-exon junctions and recover evidence for introns and transcript boundaries.

This step is done to generate transcript evidence from RNA-seq. These alignments help identify expressed regions, exon structures, and splice junctions, all of which improve gene model prediction. In a BRAKER workflow, this is one of the most important evidence-generation steps.

### `03_map_isoseq.sh`

This script maps Iso-Seq long reads to the genome, typically with a long-read aligner such as minimap2. Iso-Seq reads often span full transcripts or large parts of them, making them especially informative for transcript structure.

This step is done because long-read transcript evidence complements short-read RNA-seq. While short reads are powerful for splice-junction support, long reads are better at showing complete exon connectivity and isoform structure. Including Iso-Seq evidence can improve the biological realism of the final annotation.

### `04_qc_bams.sh`

This script performs quality control on the alignment files. That usually includes checking whether BAM files were produced correctly, whether they are sorted and indexed, and whether they are ready for downstream analysis.

This step is done because annotation pipelines are very sensitive to malformed or incomplete BAM files. A QC step helps catch problems early, such as truncated files, missing indexes, or unexpected formatting issues. It is much easier to fix these problems before running BRAKER than after an annotation job fails.

### `05_merge.sh`

This script merges RNA-seq BAM files into a combined evidence set. When multiple RNA-seq libraries, samples, or lanes are available, merging them creates a single consolidated alignment resource for annotation.

This step is done to simplify downstream evidence handling and to maximize transcript support across the genome. Instead of giving BRAKER fragmented evidence from separate runs, merging creates a broader and often more complete transcriptome signal for gene prediction.

### `06_braker.sh`

This script runs BRAKER itself. It uses the prepared genome and transcript alignment evidence to train gene prediction models and produce the final structural annotation.

This is the core step of the pipeline. BRAKER integrates ab initio prediction with RNA-supported evidence to infer gene coordinates, transcripts, coding sequences, and translated proteins. The resulting annotation files form the basis for all later downstream analyses, including functional annotation.

---

## Why this workflow is structured this way

Structural genome annotation is strongest when the evidence is prepared in a logical order. Repeat masking reduces false predictions, RNA-seq mapping provides splice evidence, Iso-Seq adds long-range transcript structure, BAM quality control ensures the evidence is usable, and BRAKER brings all these components together into a final gene set.

By separating the workflow into scripts, each stage becomes easier to:

- run independently
- debug
- rerun if needed
- document clearly
- share and reproduce across projects

---

## Final outcome

By the end of this pipeline, the repository should produce a structural annotation that describes:

- where genes are located in the genome
- how transcripts and exons are organized
- which coding sequences are predicted
- which proteins are encoded

These outputs can then be used as the starting point for a downstream **functional annotation** workflow.
