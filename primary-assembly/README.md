# Genome assembly pipeline

This repository contains a step-by-step workflow to generate a high-quality genome assembly using long-read sequencing and Hi-C scaffolding.

The goal is to go from raw sequencing data to a scaffolded genome assembly that is as complete, non-redundant, and chromosome-scale as possible.

---

## Overview

Genome assembly is a multi-step process in which raw sequencing reads are transformed into contiguous sequences and then further organized into larger scaffolds using long-range information.

This pipeline is built around the idea that a strong assembly requires both accurate contig construction and careful post-processing:

1. **Quality control of sequencing data**  
   Evaluate the raw reads before assembly to identify potential problems early.

2. **Primary assembly with HiFi reads**  
   Use high-accuracy long reads to build the initial contigs.

3. **Purging redundant haplotypic sequence**  
   Remove duplicated or alternative haplotype contigs to generate a cleaner haploid representation.

4. **Map Hi-C reads to the purged assembly**  
   Align Hi-C data to the de-redundified contigs so long-range contact information can be used for scaffolding.

5. **Scaffold with YaHS**  
   Use Hi-C contact information to order and orient contigs into larger scaffolds.

6. **Map Hi-C reads to the scaffolded assembly**  
   Realign the Hi-C data to the new scaffold structure to evaluate and refine the scaffolding.

7. **Generate contact maps with Pretext**  
   Visualize Hi-C interaction patterns to inspect scaffold quality and detect possible misassemblies.

---

## Pipeline steps

### Step 1: Quality control

The first step checks the quality of the input sequencing data. This usually includes examining read length distributions, total sequencing yield, and general sequence quality.

This step is done because assembly quality depends strongly on input data quality. Problems such as low coverage, unexpected read distributions, contamination, or poor-quality libraries can negatively affect all downstream steps. Quality control helps confirm that the data are suitable for assembly and can also guide later troubleshooting if the assembly does not behave as expected.

### Step 2: HiFi assembly with hifiasm

The second step generates the initial genome assembly using hifiasm and high-fidelity long reads. hifiasm is designed to take advantage of the accuracy and length of HiFi reads to produce highly contiguous assemblies.

This step is done because the contig-building stage is the foundation of the entire pipeline. A good primary assembly should capture as much of the genome as possible while minimizing fragmentation. HiFi-based assemblers are especially useful because they often produce contigs that are both long and accurate, making them ideal for downstream scaffolding.

### Step 3: Purging redundant sequence

The third step removes redundant contigs or haplotypic duplicates from the primary assembly. In diploid or heterozygous genomes, assemblers may retain alternative versions of the same genomic region, which can artificially inflate assembly size and complicate later analyses.

This step is done to produce a cleaner assembly that better represents a single haploid view of the genome. Purging helps reduce redundancy, improves interpretability, and makes scaffolding more reliable by avoiding duplicated representations of the same region.

### Step 4: Map Hi-C reads to the purged assembly

The fourth step aligns Hi-C read pairs to the purged contig assembly. Hi-C data provide information about physical proximity of genomic regions in the nucleus, which can be used to infer which contigs belong together and how they should be arranged.

This step is done because scaffolding requires a contact map between contigs. Before the Hi-C information can be used, the reads must be mapped to the assembly so that interaction frequencies can be calculated. Using the purged assembly at this stage helps avoid misleading contacts caused by redundant contigs.

### Step 5: Scaffold with YaHS

The fifth step uses YaHS to scaffold the purged contigs based on the mapped Hi-C reads. YaHS analyzes Hi-C contact patterns to decide how contigs should be ordered and oriented relative to one another.

This step is done to move from contigs to larger scaffold structures, ideally approaching chromosome-scale assemblies. Hi-C is especially powerful at this stage because it provides long-range linkage information that cannot be recovered from contigs alone. Scaffolding is what allows fragmented but accurate contigs to be assembled into larger genomic units.

### Step 6: Map Hi-C reads to the scaffolded assembly

The sixth step remaps the Hi-C reads to the new scaffolded assembly. After scaffolding, the structure of the reference has changed, so the Hi-C data need to be aligned again in the context of the scaffolded genome.

This step is done because evaluation of the scaffolded assembly should be based on the final scaffold structure, not only on the pre-scaffold contigs. Realigning the Hi-C data makes it possible to assess whether the scaffolding is supported by the contact signal and provides the input needed for downstream visualization and manual inspection.

### Step 7: Generate contact maps with Pretext

The final step creates Hi-C contact maps using Pretext. These maps visualize interaction frequency across the scaffolded genome and are one of the main ways to inspect scaffold quality.

This step is done because contact maps allow manual validation of the assembly. A good contact map can reveal whether scaffolds are internally consistent, whether joins are well supported, and whether there are suspicious patterns that may indicate inversions, translocations, misjoins, or other assembly issues. This visual inspection is often essential before treating the assembly as final.

---

## Why this workflow is structured this way

Genome assembly is strongest when contiguity, redundancy reduction, and long-range scaffolding are treated as separate but connected steps. The initial long-read assembly builds the sequence backbone, purging simplifies the assembly into a cleaner representation, Hi-C mapping provides long-range linkage evidence, scaffolding organizes the contigs, and contact maps allow the final structure to be checked visually.

By separating the workflow into distinct steps, the pipeline becomes easier to:

- run and troubleshoot in stages
- evaluate intermediate results
- rerun only the necessary parts when something changes
- document clearly
- reproduce across projects

---

## Final outcome

By the end of this pipeline, the repository should produce a genome assembly that is:

- built from high-accuracy long reads
- reduced in redundancy
- scaffolded using Hi-C data
- supported by contact-map inspection

This final assembly can then serve as the reference for downstream analyses such as repeat annotation, gene prediction, comparative genomics, and population genomics.
