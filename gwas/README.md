# GWAS pipeline

This repository contains a step-by-step workflow to perform a genome-wide association study (GWAS), starting from sequencing reads and ending with association statistics and Manhattan plots.

The goal is to identify genomic variants associated with a phenotype of interest by combining read mapping, SNP discovery, variant filtering and merging, association testing, and visualization.

---

## Overview

A GWAS links genetic variation to phenotypic variation across individuals. In practice, this means first generating a reliable genotype dataset, then combining it with a phenotype file, and finally testing whether any variants show a statistical association with the trait.

This pipeline is built around the idea that good GWAS results depend on both clean genotype data and a clearly formatted phenotype dataset:

1. **Map sequencing reads to the reference genome**  
   Align reads so genetic variants can be identified relative to the reference.

2. **Call SNPs**  
   Detect variable sites across individuals.

3. **Merge VCF files**  
   Combine variant information into a unified dataset for downstream analysis.

4. **Run association testing with GEMMA**  
   Test each SNP for association with the phenotype while accounting for population structure or relatedness if needed.

5. **Visualize results with Manhattan plots in R**  
   Plot association statistics across the genome to identify significant peaks.

---

## Pipeline steps

### Step 1: Read mapping

The first step aligns sequencing reads from each sample to the reference genome. This creates the positional framework needed to identify sequence differences relative to the reference.

This step is done because SNP calling requires reads to be placed accurately on the genome. Good mapping is essential for reliable variant discovery, since poor alignments can create false SNPs or obscure real ones. In a GWAS workflow, the quality of the association results depends heavily on the quality of the mapping step.

### Step 2: SNP calling

The second step identifies single nucleotide polymorphisms across the mapped samples. SNP calling compares the aligned reads to the reference genome and determines where individuals differ at specific positions.

This step is done to generate the raw genotype information used in the GWAS. Each SNP represents a candidate site that may or may not be associated with the phenotype. The goal at this stage is not yet to interpret function, but to build a consistent set of genotypes across all individuals.

### Step 3: Merge VCF files

The third step combines variant calls into a merged VCF dataset covering all samples. Depending on the workflow, this may involve merging per-sample VCFs, per-chromosome VCFs, or both.

This step is done because association analysis requires a single genotype matrix across all individuals. A merged VCF ensures that each SNP is represented consistently and that genotypes can be compared across the full sample set. This also creates a convenient point for downstream filtering, formatting, and conversion into the file types required by GWAS software.

### Step 4: Association testing with GEMMA

The fourth step runs the GWAS itself using GEMMA. GEMMA tests each SNP for association with the phenotype and can incorporate models that account for relatedness, structure, or other sources of non-independence among samples.

This step is done because simple statistical tests can be misleading when samples are related or when there is population structure. GEMMA is widely used because it provides efficient association testing and can fit linear mixed models, which are especially useful in biological datasets where individuals are not completely independent. The output of this step is a table of association statistics, usually including p-values, effect estimates, and related summary measures for each SNP.

### Step 5: Manhattan plots in R

The fifth step visualizes the GWAS results in R, most commonly with a Manhattan plot. In a Manhattan plot, each SNP is placed according to its genomic position and its statistical significance.

This step is done because visualization makes it much easier to identify association peaks and interpret the genomic distribution of signal. A Manhattan plot quickly shows whether there are strong candidate regions, whether signal is spread across the genome, or whether no variants stand out clearly. It is one of the standard outputs of a GWAS analysis.

---

## Phenotype files

A phenotype file tells the GWAS software which trait value belongs to each individual. It is one of the most important inputs in the analysis because the genotype data alone cannot be tested without a matching description of the phenotype.

### What a phenotype file contains

A phenotype file usually contains one row per individual and one or more columns representing measured traits. These traits can be:

- continuous, such as body size, wing length, or expression level
- binary, such as affected versus unaffected
- categorical, if encoded appropriately for the chosen model

In many GWAS workflows, the phenotype file must match the order of individuals in the genotype data. This is critical: if phenotype values are assigned to the wrong individuals, the association results become meaningless.

### Why phenotype files matter

The phenotype file is what defines the biological question of the GWAS. The analysis is not asking whether SNPs vary across individuals in general, but whether SNPs vary in a way that correlates with the measured trait. A carefully formatted phenotype file ensures that each genotype profile is linked to the correct biological measurement.

### Common considerations

Phenotype files often require attention to the following:

- **Sample matching**: individual names or order must correspond exactly to the genotype dataset
- **Missing data**: missing phenotypes must be coded in the format expected by the GWAS software
- **Trait type**: continuous and binary traits may require different models
- **Multiple phenotypes**: some files contain several traits, in which case the analysis may be repeated for each column

### How phenotype files are used in practice

In a typical workflow, the phenotype file is supplied to GEMMA alongside the genotype-derived input files. GEMMA then uses the phenotype values as the response variable in the association model and tests each SNP as a predictor. The quality and consistency of this phenotype file are just as important as the quality of the variant data.

---

## Why this workflow is structured this way

GWAS requires a clear progression from raw sequencing data to statistical association. Mapping creates the genomic framework, SNP calling detects variation, merging produces a unified genotype dataset, GEMMA performs the association tests, and R visualizes the results. The phenotype file links these genetic variants to the biological trait of interest.

By separating the workflow into distinct steps, the pipeline becomes easier to:

- check for errors at each stage
- rerun individual parts without restarting everything
- document and reproduce
- interpret results more clearly

---

## Final outcome

By the end of this pipeline, the repository should produce:

- a mapped read dataset
- a merged SNP dataset
- GWAS association statistics from GEMMA
- Manhattan plots showing genomic regions associated with the phenotype

Together, these outputs allow the user to identify candidate loci linked to phenotypic variation and provide a foundation for follow-up analyses and biological interpretation.
