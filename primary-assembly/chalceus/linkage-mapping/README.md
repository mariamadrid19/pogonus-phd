# Linkage mapping pipeline for *Pogonus chalceus*

This repository contains the scripts used to build a linkage map and anchor scaffolds using **Lep-MAP3** and **Lep-Anchor** for *Pogonus chalceus*.

The pipeline consists of four main steps:

1. `01_map.sh` — map sequencing reads to the reference genome
2. `02_call.sh` — generate genotype likelihoods from BAM files
3. `03_lepmap3.sh` — build linkage groups and order markers with Lep-MAP3
4. `04_lepanchor.sh` — anchor and orient scaffolds with Lep-Anchor

---

## Overview

This workflow starts from paired-end reads for mapping populations, aligns them to a fragmented reference genome, calculates posterior genotype likelihoods, constructs linkage groups with Lep-MAP3, orders markers within linkage groups, and finally uses Lep-Anchor to anchor scaffolds into chromosome-scale pseudomolecules.

This version was run assuming:

- **11 chromosomes / linkage groups**
- **achiasmatic males**, implemented with `recombination1=0` in `OrderMarkers2`

---

## Requirements

### Software

- [BWA](http://bio-bwa.sourceforge.net/)
- [SAMtools](http://www.htslib.org/)
- [Java](https://www.java.com/)
- [Lep-MAP3](https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/)
- [Lep-Anchor](https://sourceforge.net/p/lep-anchor/wiki/Home/)

### Example module versions used

- `BWA/0.7.17-GCC-10.3.0`
- `SAMtools/1.16.1-GCC-11.3.0`
- `Java/21.0.8`

---

## Input files

The pipeline expects the following main inputs:

- reference genome FASTA  
  e.g. `./genome/P_chalceus_broken.fa`
- paired-end read files  
  e.g. `./reads/SAMPLE.1.fil.fq_1.gz` and `./reads/SAMPLE.2.fil.fq_2.gz`
- sample list  
  `./samples/samples.txt`
- pedigree file for Lep-MAP3  
  `pedigree.txt`
- mapping file for `Pileup2Likelihoods`  
  `mapping.txt`
- a text file listing sorted BAM files  
  `sorted_bams.txt`

---

## Pipeline

## 1. Read mapping

Script: `01_map.sh`

This step maps paired-end reads to the reference genome using **BWA MEM**, sorts alignments, and indexes the resulting BAM files.

### Main steps

- load BWA and SAMtools
- read sample names from `samples/samples.txt`
- map reads to the reference genome
- convert SAM to BAM
- sort BAM
- index BAM

### Output

- one sorted BAM file per sample
- one BAM index per sample

Example output:

- `./bam/SAMPLE.bam`
- `./bam/SAMPLE.bam.bai`

---

## 2. Genotype likelihood calling

Script: `02_call.sh`

This step uses **samtools mpileup** and **Lep-MAP3 Pileup2Likelihoods** to calculate posterior genotype likelihoods across all mapped individuals.

### Main steps

- run `samtools mpileup`
- use `Pileup2Likelihoods` to convert pileup data into Lep-MAP3 posterior genotype format

### Output

- `post.gz`

---

## 3. Linkage mapping with Lep-MAP3

Script: `03_lepmap3.sh`

This step performs marker filtering, linkage group assignment, singleton joining, and marker ordering.

### Main steps

#### 3.1 ParentCall2

```bash
zcat post.gz | java -cp $LEPMAP/bin ParentCall2 \
    data=pedigree.txt \
    XLimit=2 \
    posteriorFile=- \
    removeNonInformative=1 | gzip > data.call.gz
