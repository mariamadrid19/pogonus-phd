library(tidyverse)
library(clusterProfiler)

setwd("~/Downloads/enrichment")

outdir <- "inversion_enrichment_results"
dir.create(outdir, showWarnings = FALSE)

# -----------------------------
# 1. Read eggNOG annotations
# -----------------------------

eggnog <- read_tsv(
  "braker_longest_isoforms.emapper.annotations",
  comment = "##",
  show_col_types = FALSE
)

names(eggnog)[1] <- "query"

# GO: TERM2GENE must be GO_term first, gene_id second
term2gene_go <- eggnog %>%
  select(gene_id = query, term = GOs) %>%
  separate_rows(term, sep = ",") %>%
  mutate(gene_id = str_remove(gene_id, "\\.t\\d+$")) %>%
  filter(!is.na(term), term != "-") %>%
  select(term, gene_id) %>%
  distinct()

# KEGG pathway enrichment, not individual KO enrichment
term2gene_kegg <- eggnog %>%
  select(gene_id = query, term = KEGG_Pathway) %>%
  separate_rows(term, sep = ",") %>%
  mutate(gene_id = str_remove(gene_id, "\\.t\\d+$")) %>%
  filter(!is.na(term), term != "-", str_starts(term, "ko")) %>%
  select(term, gene_id) %>%
  distinct()

message("Genes with GO annotations: ", n_distinct(term2gene_go$gene_id))
message("Genes with KEGG pathway annotations: ", n_distinct(term2gene_kegg$gene_id))

# -----------------------------
# 2. Read DESeq2 and inversion files
# -----------------------------

deg_all <- read_tsv("DESeq2_SW_vs_LW.tsv", show_col_types = FALSE) %>%
  filter(!is.na(padj))

inv_deg <- read_tsv("inversion_genes_ranked_total.tsv", show_col_types = FALSE) %>%
  filter(!is.na(padj))

# -----------------------------
# 3. Extract gene coordinates from BRAKER GTF
# -----------------------------

gtf <- read_tsv(
  "braker.gtf",
  col_names = FALSE,
  comment = "#",
  show_col_types = FALSE
)

gene_coords <- gtf %>%
  filter(X3 == "gene") %>%
  transmute(
    chrom = X1,
    start = X4,
    end = X5,
    gene_id = X9
  ) %>%
  distinct()

deg_all <- deg_all %>%
  left_join(gene_coords, by = "gene_id")

message("Genes in deg_all without chromosome: ", sum(is.na(deg_all$chrom)))

# -----------------------------
# 4. Enrichment function
# -----------------------------

run_enrich <- function(target_genes, background_genes, term2gene, min_genes = 10) {
  
  target_genes <- unique(na.omit(target_genes))
  background_genes <- unique(na.omit(background_genes))
  
  if (length(target_genes) < min_genes) return(NULL)
  
  enricher(
    gene = target_genes,
    universe = background_genes,
    TERM2GENE = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
}

# -----------------------------
# 5. Build chromosome/ecotype gene sets
# -----------------------------

gene_sets <- inv_deg %>%
  mutate(direction = if_else(log2FoldChange > 0, "SW_up", "LW_up")) %>%
  group_by(chrom, direction) %>%
  summarise(
    genes = list(unique(gene_id)),
    n_genes = n_distinct(gene_id),
    .groups = "drop"
  )

write_tsv(
  gene_sets %>% select(chrom, direction, n_genes),
  file.path(outdir, "gene_counts_per_chr_direction.tsv")
)

print(gene_sets %>% select(chrom, direction, n_genes))

# -----------------------------
# 6. Run GO and KEGG enrichment
# -----------------------------

results <- list()

for (i in seq_len(nrow(gene_sets))) {
  
  chr <- gene_sets$chrom[i]
  direction <- gene_sets$direction[i]
  target <- gene_sets$genes[[i]]
  
  background <- deg_all %>%
    filter(chrom == chr) %>%
    pull(gene_id)
  
  message(chr, " ", direction, ": target = ", length(target),
          "; background = ", length(background))
  
  results[[paste(chr, direction, "GO", sep = "_")]] <-
    run_enrich(target, background, term2gene_go)
  
  results[[paste(chr, direction, "KEGG", sep = "_")]] <-
    run_enrich(target, background, term2gene_kegg)
}

# -----------------------------
# 7. Combine and save all results
# -----------------------------

all_results <- bind_rows(
  lapply(names(results), function(test_name) {
    
    res <- results[[test_name]]
    if (is.null(res)) return(NULL)
    
    df <- as.data.frame(res)
    if (nrow(df) == 0) return(NULL)
    
    df %>%
      mutate(test = test_name, .before = 1)
  })
)

write_tsv(
  all_results,
  file.path(outdir, "all_enrichment_results.tsv")
)

# Save significant only
sig_results <- all_results %>%
  filter(p.adjust < 0.05)

write_tsv(
  sig_results,
  file.path(outdir, "significant_enrichment_results_FDR005.tsv")
)

# Save one file per test
for (test_name in unique(all_results$test)) {
  all_results %>%
    filter(test == test_name) %>%
    arrange(p.adjust) %>%
    write_tsv(file.path(outdir, paste0(test_name, ".tsv")))
}

# -----------------------------
# 8. Quick summaries
# -----------------------------

top_terms <- all_results %>%
  arrange(p.adjust) %>%
  group_by(test) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  select(test, ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue)

write_tsv(
  top_terms,
  file.path(outdir, "top10_terms_per_test.tsv")
)

print(top_terms, n = 100)
