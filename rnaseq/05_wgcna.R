
# -----------------------------
# WGCNA analysis
# -----------------------------
library(WGCNA)
library(data.table)


# -----------------------------
# Paths
# -----------------------------
vst_file <- "deseq2/vst_counts.tsv"
samples_file <- "samples.tsv"
annot_file <- "braker.minimal.annotation_map.tsv"
outdir <- "wgcna"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(outdir, "module_gene_lists"), showWarnings = FALSE, recursive = TRUE)


# -----------------------------
# Read data
# -----------------------------
vst <- fread(vst_file, data.table = FALSE, check.names = FALSE)
samples <- fread(samples_file, data.table = FALSE)
annot <- fread(annot_file, data.table = FALSE)

rownames(vst) <- vst$gene_id
vst$gene_id <- NULL

# Keep only samples in metadata, in metadata order
vst <- vst[, samples$sample, drop = FALSE]

# WGCNA expects rows = samples, columns = genes
datExpr0 <- as.data.frame(t(vst))

# -----------------------------
# Trait data
# -----------------------------
samples$condition <- factor(samples$condition, levels = c("LW", "SW"))
samples$condition_num <- ifelse(samples$condition == "SW", 1, 0)
rownames(samples) <- samples$sample

traitData <- samples[, c("condition_num"), drop = FALSE]
colnames(traitData) <- c("SW_vs_LW")

# -----------------------------
# Check genes/samples
# -----------------------------
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    message("Removing genes: ", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", "))
  }
  if (sum(!gsg$goodSamples) > 0) {
    message("Removing samples: ", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", "))
  }
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]
  traitData <- traitData[rownames(datExpr0), , drop = FALSE]
}

# -----------------------------
# Filter low-variance genes
# -----------------------------
gene_vars <- apply(datExpr0, 2, var, na.rm = TRUE)

# Keep top 5000 variable genes, or all if fewer
n_keep <- min(5000, length(gene_vars))
keep_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_keep]

datExpr <- datExpr0[, keep_genes, drop = FALSE]

write.table(
  data.frame(gene_id = keep_genes, variance = gene_vars[keep_genes]),
  file = file.path(outdir, "genes_used_for_wgcna.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Sample clustering
# -----------------------------
sampleTree <- hclust(dist(datExpr), method = "average")

pdf(file.path(outdir, "01_sample_clustering.pdf"), width = 8, height = 6)
par(cex = 0.8)
plot(sampleTree, main = "Sample clustering", xlab = "", sub = "")
dev.off()

# -----------------------------
# Soft-threshold selection
# -----------------------------
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

pdf(file.path(outdir, "02_soft_threshold.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))

cex1 <- 0.9
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1
)
abline(h = 0.8, col = "red")

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1
)
dev.off()

# Choose power: first >= 0.8 if possible, otherwise 6
fit_df <- sft$fitIndices
signed_r2 <- -sign(fit_df[, 3]) * fit_df[, 2]
candidate_powers <- fit_df[signed_r2 >= 0.8, 1]

if (length(candidate_powers) > 0) {
  softPower <- min(candidate_powers)
} else {
  softPower <- 6
}

writeLines(
  paste("Chosen soft-threshold power:", softPower),
  con = file.path(outdir, "chosen_soft_power.txt")
)

# -----------------------------
# Network construction and module detection
# -----------------------------
net <- blockwiseModules(
  datExpr,
  power = softPower,
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = TRUE,
  saveTOMs = FALSE,
  verbose = 3
)

moduleColors <- net$colors
moduleLabels <- moduleColors
MEs <- net$MEs

# Plot gene dendrogram and module colors
pdf(file.path(outdir, "03_gene_dendrogram_modules.pdf"), width = 12, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
dev.off()

# -----------------------------
# Module-trait relationships
# -----------------------------
MEs0 <- moduleEigengenes(datExpr, colors = moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

moduleTraitCor <- cor(MEs, traitData, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

textMatrix <- paste(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 2),
  ")",
  sep = ""
)

pdf(file.path(outdir, "04_module_trait_heatmap.pdf"), width = 7, height = 8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(traitData),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module-trait relationships"
)
dev.off()

# Save module-trait table
module_trait_df <- data.frame(
  module = rownames(moduleTraitCor),
  trait = colnames(traitData)[1],
  correlation = moduleTraitCor[, 1],
  pvalue = moduleTraitPvalue[, 1],
  stringsAsFactors = FALSE
)
module_trait_df <- module_trait_df[order(module_trait_df$pvalue), ]

write.table(
  module_trait_df,
  file = file.path(outdir, "module_trait_relationships.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Gene-module membership table
# -----------------------------
geneInfo <- data.frame(
  gene_id = colnames(datExpr),
  module = moduleColors,
  stringsAsFactors = FALSE
)

# Module membership (kME)
for (module in substring(names(MEs), 3)) {
  column <- match(paste0("ME", module), colnames(MEs))
  geneInfo[[paste0("kME_", module)]] <- cor(datExpr, MEs[, column], use = "p")
}

# Gene significance for trait
geneTraitSignificance <- cor(datExpr, traitData$SW_vs_LW, use = "p")
geneTraitPvalue <- corPvalueStudent(geneTraitSignificance, nrow(datExpr))
geneInfo$GS_SW_vs_LW <- as.numeric(geneTraitSignificance)
geneInfo$GS_pvalue <- as.numeric(geneTraitPvalue)

# Join annotation
geneInfo <- merge(geneInfo, annot, by = "gene_id", all.x = TRUE, sort = FALSE)

write.table(
  geneInfo,
  file = file.path(outdir, "gene_module_membership.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# -----------------------------
# Write per-module gene lists
# -----------------------------
for (mod in sort(unique(moduleColors))) {
  mod_df <- subset(geneInfo, module == mod)
  mod_df <- mod_df[order(-abs(mod_df$GS_SW_vs_LW)), ]
  write.table(
    mod_df,
    file = file.path(outdir, "module_gene_lists", paste0("module_", mod, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# -----------------------------
# Plot eigengenes for top modules
# -----------------------------
top_modules <- head(module_trait_df$module, 6)

for (me_name in top_modules) {
  if (!me_name %in% colnames(MEs)) next
  
  plot_df <- data.frame(
    sample = rownames(MEs),
    eigengene = MEs[, me_name],
    condition = samples[rownames(MEs), "condition"]
  )
  
  p <- ggplot(plot_df, aes(x = condition, y = eigengene)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0) +
    theme_bw(base_size = 12) +
    labs(
      title = paste("Module eigengene:", me_name),
      x = "Condition",
      y = "Eigengene"
    )
  
  ggsave(
    filename = file.path(outdir, paste0("eigengene_", me_name, ".pdf")),
    plot = p,
    width = 5,
    height = 4
  )
}

# -----------------------------
# Heatmap of top module eigengenes
# -----------------------------
top_me_mat <- MEs[, top_modules, drop = FALSE]
annotation_col <- data.frame(condition = samples[rownames(top_me_mat), "condition"])
rownames(annotation_col) <- rownames(top_me_mat)

pdf(file.path(outdir, "05_top_module_eigengenes_heatmap.pdf"), width = 7, height = 5)
pheatmap(
  t(top_me_mat),
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Top module eigengenes"
)
dev.off()

# -----------------------------
# Summary file
# -----------------------------
sink(file.path(outdir, "summary.txt"))
cat("WGCNA summary\n")
cat("====================\n")
cat("Samples used:", nrow(datExpr), "\n")
cat("Genes used:", ncol(datExpr), "\n")
cat("Chosen soft-threshold power:", softPower, "\n\n")

cat("Module sizes:\n")
print(sort(table(moduleColors), decreasing = TRUE))

cat("\nTop module-trait relationships:\n")
print(head(module_trait_df, 10))
sink()

message("WGCNA finished. Results written to: ", outdir)
