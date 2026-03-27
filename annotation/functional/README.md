# Functional annotation pipeline for BRAKER gene predictions

This repository contains a step-by-step functional annotation workflow for a genome annotated with **BRAKER**. The goal is to start from BRAKER structural gene predictions and build a reproducible pipeline that assigns functional information to predicted genes.

This workflow is designed around the following BRAKER output files:

- `braker.aa` — predicted protein sequences
- `braker.codingseq` — coding nucleotide sequences (CDS)
- `braker.gtf` — gene model coordinates and transcript structure

In this project, the files are stored here:

- BRAKER input directory: `/scratch/leuven/357/vsc35707/annotation/braker`
- Functional annotation project directory: `/scratch/leuven/357/vsc35707/annotation/func-annotation`

The pipeline is built so that everything related to functional annotation lives under `func-annotation/`.

---

## Why this pipeline starts from BRAKER output

BRAKER performs **structural annotation**: it predicts where genes, transcripts, exons, and CDS regions are located in the genome. However, BRAKER does **not** tell you in a biologically meaningful way what those genes do.

Functional annotation is the downstream process that tries to answer questions like:

- What type of protein does this gene encode?
- Does it have a known conserved domain?
- Does it resemble a known gene from another species?
- Can GO terms, pathways, or orthology assignments be transferred?

The general strategy used here is:

1. clean and standardize the BRAKER outputs
2. keep one representative isoform per gene
3. assess overall gene set quality
4. assign orthology-based function
5. assign domain-based function
6. optionally add curated similarity-based hits
7. merge all evidence into one final annotation table

---

## Recommended project structure

```text
/scratch/leuven/357/vsc35707/annotation/func-annotation
├── config
├── input
│   └── braker
├── logs
├── results
│   ├── 01_longest_isoforms
│   ├── 02_busco
│   ├── 03_eggnog
│   ├── 04_interproscan
│   ├── 05_diamond
│   └── 06_final_annotation
├── scripts
├── tmp
└── databases
    └── eggnog
```

Suggested logic for these directories:

- `input/braker/` contains symlinks to the BRAKER output files
- `scripts/` contains one script per analysis step
- `logs/` stores SLURM stdout/stderr logs
- `results/` stores outputs from each major step
- `tmp/` stores temporary run files
- `databases/` stores local annotation databases such as eggNOG

---

## Step 0: Link BRAKER output files into the project

Instead of moving the original BRAKER files, this workflow uses symlinks. This keeps the original structural annotation directory unchanged while making all downstream paths consistent.

Example:

```bash
BRAKER=/scratch/leuven/357/vsc35707/annotation/braker
PROJECT=/scratch/leuven/357/vsc35707/annotation/func-annotation

ln -s $BRAKER/braker.aa        $PROJECT/input/braker/braker.aa
ln -s $BRAKER/braker.codingseq $PROJECT/input/braker/braker.codingseq
ln -s $BRAKER/braker.gtf       $PROJECT/input/braker/braker.gtf
```

---

## Step 1: Extract the longest isoform per gene

### Why this step matters

BRAKER often predicts multiple transcript isoforms per gene. Running functional annotation on all isoforms can create redundancy and make downstream summaries harder to interpret. For most whole-genome annotation workflows, it is standard to keep **one representative protein per gene**, usually the **longest isoform**.

This gives a cleaner, less redundant protein set for tools like:

- BUSCO
n- eggNOG-mapper
- InterProScan
- DIAMOND

### Script

- `scripts/extract_longest_isoforms.py`

### Inputs

- `input/braker/braker.gtf`
- `input/braker/braker.aa`
- `input/braker/braker.codingseq`

### Outputs

Written to `results/01_longest_isoforms/`:

- `braker.longest_isoforms.aa`
- `braker.longest_isoforms.codingseq`
- `braker.longest_isoforms.tsv`

### What the script does

The script:

1. parses the GTF file to map `transcript_id` to `gene_id`
2. reads the protein FASTA
3. calculates the protein length of each transcript
4. selects the longest transcript per gene
5. writes filtered FASTA files for proteins and CDS
6. writes a summary table recording which transcript was kept for each gene

### Result obtained in this run

- original protein set: `17,110` proteins
- representative protein set after filtering: `13,531` proteins

This indicates that a substantial number of alternative isoforms were removed, which is expected and desirable for downstream annotation.

### Notes

- This step assumes the FASTA headers use transcript IDs compatible with the `transcript_id` values in the GTF.
- In case of equal-length isoforms, the script uses a deterministic tie-break rule.

---

## Step 2: Assess gene set completeness with BUSCO

### Why this step matters

Before spending time on functional annotation, it is important to verify that the predicted proteome is reasonably complete. BUSCO evaluates the presence of highly conserved single-copy orthologs expected for a taxonomic lineage.

This provides a quick quality check on whether the predicted gene set is biologically plausible and sufficiently complete.

### Script

- `scripts/run_busco_longest_isoforms.sh`

### Input

- `results/01_longest_isoforms/braker.longest_isoforms.aa`

### Output directory

- `results/02_busco/`

### Lineage used

- `insecta_odb12`

This is appropriate as a first-pass lineage for a beetle proteome.

### Result obtained in this run

BUSCO summary:

- `C:96.1%`
- `S:94.4%`
- `D:1.7%`
- `F:0.9%`
- `M:3.0%`
- `n:3114`

### Interpretation

This is a strong result:

- high completeness suggests good gene-space recovery
- high single-copy proportion suggests low redundancy after isoform filtering
- low duplication suggests the representative set is clean
- low fragmented/missing fraction suggests the proteome is suitable for downstream functional annotation

Because this BUSCO result is strong, it is reasonable to proceed with annotation.

---

## Step 3: Annotate proteins with eggNOG-mapper

### Why this step matters

eggNOG-mapper provides **orthology-based functional annotation**. It is one of the most useful first-pass annotation tools for a new genome because it can assign:

- orthologs
- gene descriptions
- GO terms
- functional categories
- KEGG-related annotations
- broader evolutionary context

This step is especially useful for summarizing the likely identity of predicted proteins.

### Script

- `scripts/03_run_eggnog.sh`

### Input

- `results/01_longest_isoforms/braker.longest_isoforms.aa`

### Output directory

- `results/03_eggnog/`

### Database location

- `/scratch/leuven/357/vsc35707/annotation/func-annotation/databases/eggnog`

### Software environment

The run assumes:

- `emapper.py`
- `diamond`
- conda environment: `eggnog`

### Important database note

The standard `download_eggnog_data.py` script in this setup pointed to an outdated eggNOG domain and returned HTTP 404 errors. The databases were therefore downloaded manually.

For DIAMOND-based annotation, the required files are:

- `eggnog.db`
- `eggnog.taxa.db`
- `eggnog_proteins.dmnd`

A beetle-specific or Coleoptera-specific download is **not required** for the normal DIAMOND workflow. Taxonomic restriction is instead handled during annotation, for example with:

```bash
--tax_scope eukaryota
```

### What the script does

The script:

1. activates the `eggnog` conda environment
2. checks that input and database files exist
3. runs `emapper.py` in protein mode with DIAMOND
4. writes all outputs into `results/03_eggnog/`
5. removes temporary files after completion

### Main output expected

The most important file will usually be:

- `braker_longest_isoforms.emapper.annotations`

This table will later be merged with other evidence.

---

## Step 4: Annotate conserved domains with InterProScan

### Why this step matters

Sequence similarity alone is not always enough to give a trustworthy annotation. Domain-based annotation provides an independent line of evidence.

InterProScan integrates many protein signature databases and can identify:

- protein domains
- conserved sites
- repeats
- family assignments
- GO terms linked to domain evidence

This is especially valuable when:

- eggNOG hits are weak
- proteins are lineage-specific
- a protein is clearly real but lacks a strong ortholog assignment

### Planned script

- `scripts/04_run_interproscan.sh`

### Expected input

- `results/01_longest_isoforms/braker.longest_isoforms.aa`

### Expected output directory

- `results/04_interproscan/`

### Typical useful outputs

- InterProScan TSV
- InterProScan GFF3
- GO term assignments
- pathway-related output if enabled

### Role in the pipeline

InterProScan adds structural and biochemical context that complements eggNOG orthology calls.

---

## Step 5: Similarity search against a curated protein database

### Why this step matters

A curated similarity search provides an additional evidence layer for naming proteins. The preferred target is typically a reviewed protein database such as Swiss-Prot.

This step helps with:

- assigning conservative protein names
- spotting likely contaminants or odd predictions
- manually checking genes of special biological interest

### Planned script

- `scripts/05_run_diamond_swissprot.sh`

### Expected input

- `results/01_longest_isoforms/braker.longest_isoforms.aa`

### Expected output directory

- `results/05_diamond/`

### Role in the pipeline

This step is best used as a naming and manual curation layer rather than the only annotation source.

---

## Step 6: Merge evidence into a final annotation table

### Why this step matters

No single annotation method should be trusted on its own. The strongest annotations usually come from combining multiple evidence types.

The final annotation table should ideally include columns such as:

- gene ID
- transcript ID
- protein length
- selected representative isoform
- eggNOG preferred name
- eggNOG description
- GO terms
- KEGG terms
- InterPro accession(s)
- InterPro description(s)
- domain/family information
- best curated DIAMOND hit
- percent identity
- e-value
- confidence tier

### Planned script

- `scripts/06_merge_annotations.py`

### Expected output directory

- `results/06_final_annotation/`

### Suggested confidence logic

A practical confidence ranking can be:

- **high confidence**: strong curated similarity hit and compatible InterPro domains
- **medium confidence**: eggNOG orthology and/or clear InterPro domain evidence
- **low confidence**: weak or generic similarity only
- **uncharacterized**: no convincing support

### Role in the pipeline

This merged table is the main deliverable for downstream biology.

---

## Step 7: Add functional labels back onto the genome annotation

### Why this step matters

The final goal is not only to have a spreadsheet of annotations, but also to be able to associate functional labels with the coordinates in the BRAKER annotation.

This can be done by joining the annotation table back to the original GTF/GFF identifiers and adding attributes such as:

- product name
- note/description
- GO terms
- InterPro IDs
- external database references
