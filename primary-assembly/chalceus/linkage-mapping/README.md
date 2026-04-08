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
- paired-end read files  
- sample list  
- pedigree file for Lep-MAP3  
- mapping file for `Pileup2Likelihoods`  
- a text file listing sorted BAM files  

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
zcat post.gz | java -cp $LEPMAP/bin ParentCall2 data=pedigree.txt XLimit=2 posteriorFile=- removeNonInformative=1 | gzip > data.call.gz
```

#### 3.2 Filtering2

```bash
zcat data.call.gz | java -cp $LEPMAP/bin Filtering2 data=- dataTolerance=0.01 | gzip > data_f_t01.call.gz
```

#### 3.3 Extract SNP names

```bash
zcat data_f_t01.call.gz | awk 'NR>=7' | cut -f 1,2 > snps.txt
```

#### 3.4 Test multiple `lodLimit` values with `SeparateChromosomes2`

```bash
for lod in 4 5 6 7 8 9 10 11 12 15; do
    zcat data_f_t01.call.gz | java -cp $LEPMAP/bin SeparateChromosomes2 data=- lodLimit=$lod > map${lod}.txt
done
```

This was used to evaluate linkage group fragmentation and identify a suitable `lodLimit`.

#### 3.5 Test multiple `lodLimit` values with `JoinSingles2All`

```bash
for lod in 4 5 6 7 8; do
    zcat data_f_t01.call.gz | java -cp $LEPMAP/bin JoinSingles2All map=map12.txt data=- lodLimit=$lod iterate=2 > map12_js${lod}.txt
done
```

#### 3.6 Order markers

Markers were ordered for each chromosome with the achiasmatic male model:

```bash
for chr in {1..11}; do
    zcat data_f_t01.call.gz | java -cp $LEPMAP/bin OrderMarkers2 map=map12.txt data=- recombination1=0 chromosome=$chr > map12_chr${chr}_mrecom0.txt
done
```

### Main output files

- `data.call.gz`
- `data_f_t01.call.gz`
- `snps.txt`
- `map*.txt`
- `map12_js*.txt`
- `map12_chr*_mrecom0.txt`

---

## 4. Scaffold anchoring with Lep-Anchor

Script: `04_lepanchor.sh`

This step uses linkage map information to clean map assignments, convert them to BED intervals, and place/orient contigs using Lep-Anchor.

### User-defined variables

```bash
MAP_PREFIX="map12"
NCHR=11
USE_MRECOM0=1
```

These variables define:

- the linkage map prefix to use
- the number of chromosomes
- whether to use the `_mrecom0` marker order files

### Main steps

#### 4.1 Count markers per linkage group

```bash
cut -f1 "$MAP_FILE" | sort -n | uniq -c
```

#### 4.2 Prepare CleanMap input

```bash
paste snps.txt "$MAP_FILE" | awk '(NR>1)' > "$CLEANMAP_INPUT"
sort -V -k1,1 -k2,2n "$CLEANMAP_INPUT" > "$CLEANMAP_SORTED"
```

#### 4.3 Run `CleanMap`

```bash
java -cp "$LEPANCHOR/bin/" CleanMap map="$CLEANMAP_SORTED" > "$MAP_CLEAN"
```

#### 4.4 Generate genome size file

```bash
awk -f "$LEPANCHOR/contigLength.awk" "$GENOME" > "$GENOME.sizes"
```

#### 4.5 Convert cleaned map to BED

```bash
java -cp "$LEPANCHOR/bin/" Map2Bed map="$MAP_CLEAN" contigLength="$GENOME.sizes" > "$MAP_BED"
```

#### 4.6 Prepare `PlaceAndOrientContigs` input

```bash
for X in $(seq 1 "$NCHR"); do
  awk -vn="$X" '
    (NR==FNR){map[NR-1]=$0}
    (NR!=FNR && /^[^#]/){print map[$1], n, $2, $3}
  ' snps.txt "${MAP_PREFIX}_chr${X}${ORDER_SUFFIX}.txt" > "${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.input"
done
```

#### 4.7 Run `PlaceAndOrientContigs`

```bash
for X in $(seq 1 "$NCHR"); do
  java -cp "$LEPANCHOR/bin/" PlaceAndOrientContigs map="${MAP_PREFIX}_chr${X}${MINPUT_SUFFIX}.input" bed="$MAP_BED" chromosome="$X" noIntervals=1 > "${MAP_PREFIX}_chr${X}.la" 2> "${MAP_PREFIX}_chr${X}.la.err"
done
```

#### 4.8 Generate AGP files

```bash
for X in $(seq 1 "$NCHR"); do
  awk -vlg="$X" -f "$LEPANCHOR/makeagp_full2.awk" "${MAP_PREFIX}_chr${X}.la" > "${MAP_PREFIX}_chr${X}.agp"
done
```

### Main output files

- `*.clean`
- `*.bed`
- `*_chr*.input`
- `*_chr*.la`
- `*_chr*.la.err`
- `*_chr*.agp`

---

## Final FASTA generation

After generating AGP files, a chromosome-scale FASTA can be created with:

```bash
awk -f /data/leuven/357/vsc35707/LepAnchor/makefasta.awk /scratch/leuven/357/vsc35707/linkage-mapping/genome/Pchalceus_LW.fa chr*.agp > final.fasta
```

If desired, unanchored scaffolds can then be appended to produce a final assembly containing both chromosomes and extra scaffolds.

---

## Notes

- The best `lodLimit` for `SeparateChromosomes2` was assessed empirically by comparing marker distributions across linkage groups.
- `recombination1=0` was used because male carabid beetles are assumed to be **achiasmatic**.
- `JoinSingles2All` was explored with multiple `lodLimit` values to balance marker rescue and over-assignment.
- `CleanMap` input must be sorted by scaffold and position before running Lep-Anchor.
