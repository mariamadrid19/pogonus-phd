library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)

counts_file <- "gene_counts.txt"
samples_file <- "samples.tsv"
outdir <- "deseq2"


# Read counts
counts <- read.delim(counts_file, comment.char = "#", check.names = FALSE)

# Keep gene IDs and count columns only
count_matrix <- counts[, c("Geneid", grep("\\.filtered\\.sorted\\.bam$", colnames(counts), value = TRUE))]

rownames(count_matrix) <- count_matrix$Geneid
count_matrix <- count_matrix[, -1]

# Clean column names to match sample sheet
colnames(count_matrix) <- basename(colnames(count_matrix))
colnames(count_matrix) <- sub("\\.filtered\\.sorted\\.bam$", "", colnames(count_matrix))

# Read metadata
coldata <- read.delim(samples_file, stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample

# Reorder metadata to match count matrix
coldata <- coldata[colnames(count_matrix), , drop = FALSE]
coldata$condition <- factor(coldata$condition, levels = c("LW", "SW"))

# Build DESeq object
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_matrix)),
  colData = coldata,
  design = ~ condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Run DESeq2
dds <- DESeq(dds)

# Main contrast: SW vs LW
res <- results(dds, contrast = c("condition", "SW", "LW"))
res <- lfcShrink(dds, coef = "condition_SW_vs_LW", type = "apeglm")

# Order by adjusted p-value
res_ordered <- res[order(res$padj), ]
res_df <- as.data.frame(res_ordered)
res_df$gene_id <- rownames(res_df)

write.table(
  res_df,
  file = file.path(outdir, "DESeq2_SW_vs_LW.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
norm_counts <- as.data.frame(norm_counts)
norm_counts$gene_id <- rownames(norm_counts)

write.table(
  norm_counts,
  file = file.path(outdir, "normalized_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# VST for PCA
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
vst_df <- as.data.frame(vst_mat)
vst_df$gene_id <- rownames(vst_df)

write.table(
  vst_df,
  file = file.path(outdir, "vst_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# PCA plot
pdf(file.path(outdir, "PCA.pdf"))
plotPCA(vsd, intgroup = "condition")
dev.off()

# MA plot
pdf(file.path(outdir, "MAplot.pdf"))
plotMA(res, ylim = c(-5, 5))
dev.off()

# Summary
sink(file.path(outdir, "summary.txt"))
cat("Samples:\n")
print(coldata)
cat("\nDESeq2 results summary:\n")
print(summary(res))
sink()



annot_file <- "braker.minimal.annotation_map.tsv"
deg_file <- "deseq2/DESeq2_SW_vs_LW.tsv"

annot <- read.delim(annot_file, stringsAsFactors = FALSE)
deg <- read.delim(deg_file, stringsAsFactors = FALSE)

deg_annot <- merge(deg, annot, by = "gene_id", all.x = TRUE, sort = FALSE)

write.table(
  deg_annot,
  file = "DESeq2_SW_vs_LW.annotated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


sig <- subset(deg_annot, !is.na(padj) & padj < 0.1)
up  <- subset(sig, log2FoldChange > 0)
down <- subset(sig, log2FoldChange < 0)

#write.table(sig,file = "deseq2/DESeq2_SW_vs_LW.annotated.sig.tsv",sep = "\t", quote = FALSE, row.names = FALSE)

#write.table(up,file = "deseq2/DESeq2_SW_vs_LW.annotated.up.tsv",sep = "\t", quote = FALSE, row.names = FALSE)

#write.table(down,file = "deseq2/DESeq2_SW_vs_LW.annotated.down.tsv",sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# Paths
# -----------------------------
deg_file   <- "DESeq2_SW_vs_LW.annotated.tsv"
vst_file   <- "deseq2/vst_counts.tsv"
samples_file <- "samples.tsv"
outdir     <- "deseq2/plots"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Read data
# -----------------------------
deg <- read.delim(deg_file, stringsAsFactors = FALSE, check.names = FALSE)
samples <- read.delim(samples_file, stringsAsFactors = FALSE)
rownames(samples) <- samples$sample

vst <- read.delim(vst_file, stringsAsFactors = FALSE, check.names = FALSE)
rownames(vst) <- vst$gene_id
vst <- vst[, setdiff(colnames(vst), "gene_id"), drop = FALSE]

# Make sure sample order matches
vst <- vst[, samples$sample, drop = FALSE]

# -----------------------------
# Volcano plot
# -----------------------------
deg$neglog10_padj <- -log10(deg$padj)
deg$neglog10_padj[is.infinite(deg$neglog10_padj)] <- NA

deg$status <- "Not significant"
deg$status[!is.na(deg$padj) & deg$padj < 0.1 & deg$log2FoldChange > 0] <- "Up in SW"
deg$status[!is.na(deg$padj) & deg$padj < 0.1 & deg$log2FoldChange < 0] <- "Up in LW"

# Label a few top genes
deg_label <- deg[!is.na(deg$padj), ]
deg_label <- deg_label[order(deg_label$padj), ]
deg_label <- head(deg_label, 20)

label_col <- if ("gene_name" %in% colnames(deg)) "gene_name" else "gene_id"
deg_label$label <- ifelse(
  !is.na(deg_label[[label_col]]) & deg_label[[label_col]] != "",
  deg_label[[label_col]],
  deg_label$gene_id
)

pdf(file.path(outdir, "volcano_SW_vs_LW.pdf"), width = 8, height = 7)
p <- ggplot(deg, aes(x = log2FoldChange, y = neglog10_padj, color = status)) +
  geom_point(alpha = 0.7, size = 1.2, na.rm = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
  labs(
    title = "Volcano plot: SW vs LW",
    x = "log2 fold change",
    y = expression(-log[10](adjusted~pvalue))
  ) +
  theme_bw(base_size = 12)

if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(
    data = deg_label,
    aes(label = label),
    size = 3,
    max.overlaps = 50,
    show.legend = FALSE
  )
}

print(p)
dev.off()

# -----------------------------
# Heatmap: top significant genes
# -----------------------------
sig <- subset(deg, !is.na(padj) & padj < 0.1)

# Pick top 50 by adjusted p-value
sig <- sig[order(sig$padj), ]
top_n <- min(50, nrow(sig))
top_sig <- head(sig, top_n)

heat_ids <- top_sig$gene_id
heat_mat <- as.matrix(vst[heat_ids, , drop = FALSE])

# Replace row names with gene_name when available and non-empty
if ("gene_name" %in% colnames(top_sig)) {
  row_labels <- ifelse(
    !is.na(top_sig$gene_name) & top_sig$gene_name != "",
    paste(top_sig$gene_name, top_sig$gene_id, sep = " | "),
    top_sig$gene_id
  )
} else {
  row_labels <- top_sig$gene_id
}
rownames(heat_mat) <- row_labels

# Z-score by gene
heat_mat_scaled <- t(scale(t(heat_mat)))

annotation_col <- data.frame(condition = samples$condition)
rownames(annotation_col) <- samples$sample

pdf(file.path(outdir, "heatmap_top50_sig_genes_SW_vs_LW.pdf"), width = 8, height = 10)
pheatmap(
  heat_mat_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  main = "Top 50 significant genes"
)
dev.off()

# -----------------------------
# Heatmap: top variable genes among significant ones
# -----------------------------
if (nrow(sig) > 2) {
  sig_ids <- sig$gene_id
  sig_mat <- as.matrix(vst[sig_ids, , drop = FALSE])
  vars <- apply(sig_mat, 1, var, na.rm = TRUE)
  vars <- sort(vars, decreasing = TRUE)
  
  top_var_n <- min(100, length(vars))
  top_var_ids <- names(vars)[1:top_var_n]
  
  var_mat <- as.matrix(vst[top_var_ids, , drop = FALSE])
  var_mat_scaled <- t(scale(t(var_mat)))
  
  sig2 <- sig[match(top_var_ids, sig$gene_id), , drop = FALSE]
  if ("gene_name" %in% colnames(sig2)) {
    row_labels2 <- ifelse(
      !is.na(sig2$gene_name) & sig2$gene_name != "",
      paste(sig2$gene_name, sig2$gene_id, sep = " | "),
      sig2$gene_id
    )
  } else {
    row_labels2 <- sig2$gene_id
  }
  rownames(var_mat_scaled) <- row_labels2
  
  pdf(file.path(outdir, "heatmap_top100_variable_sig_genes_SW_vs_LW.pdf"), width = 8, height = 12)
  pheatmap(
    var_mat_scaled,
    annotation_col = annotation_col,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 5,
    main = "Top 100 variable significant genes"
  )
  dev.off()
}

# -----------------------------
# Save top genes table used in plots
# -----------------------------
write.table(
  top_sig,
  file = file.path(outdir, "top50_sig_genes_used_for_heatmap.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Plots written to:", outdir, "\n")
