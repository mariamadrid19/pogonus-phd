library(data.table)
library("SNPRelate")
library("ggplot2")
library(car)
library(nlme)
library(emmeans)
library(effects)

setwd("/scratch/leuven/357/vsc35707/GWAS/NP25")
chr3_phased_inversion.vcf.gz

# Define input and output files
vcfFile <- "chr3_phased_inversion.vcf.gz"
gdsFile <- "chr3_phased_inversion.gds"

# Convert VCF to GDS format
snpgdsVCF2GDS(vcfFile, gdsFile, method = "biallelic.only")

# Open GDS file for analysis
genofile <- snpgdsOpen(gdsFile)

# Perform PCA analysis on the genotype data
ccm_pca <- snpgdsPCA(genofile, autosome.only = TRUE)

# Extract sample names and clean them by removing extra suffixes (whatever you have in your sample names, you maybe cleaned them already)
sNames <- ccm_pca$sample.id
sNames <- sub(".BarSW.filtered.sorted.bam", "", sNames)

popInfo <- data.frame(
  V1 = c("NP25")
)

# Create a data frame for plotting PCA results
plot_data <- as.data.frame(ccm_pca$eigenvect)
colnames(plot_data) <- paste0("PC", 1:ncol(plot_data))  # Rename columns for clarity

# divide in three clusters based on PC1
plot_data$karyotype <- cut(plot_data$PC1, breaks = 3, labels = c("type1", "type2", "type3"))
plot_data$IID <- ccm_pca$sample.id

# PCA Plot Colored by Population
ggplot(plot_data, aes(x = PC1, y = PC2, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Population", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# PC2 and PC3
ggplot(plot_data, aes(x = PC2, y = PC3, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Population", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# IID lists of different karyotypes
iid_lists <- split(plot_data$IID, plot_data$karyotype)
iid_lists$type1
iid_lists$type2
iid_lists$type3

write.table(iid_lists$type1, "homozygotesLW.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(iid_lists$type2, "heterozygotes", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(iid_lists$type3, "homozygotesSW.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


get_combinations <- function(iid_chr2, iid_chr3) {
  all_combos <- list()
  
  for (t2 in names(iid_chr2)) {
    for (t3 in names(iid_chr3)) {
      inds_chr2 <- iid_chr2[[t2]]
      inds_chr3 <- iid_chr3[[t3]]
      
      # Get intersection (individuals present in both)
      common_inds <- intersect(inds_chr2, inds_chr3)
      
      combo_name <- paste(t2, t3, sep = "_")
      all_combos[[combo_name]] <- common_inds
    }
  }
  
  return(all_combos)
}

combo_table <- get_combinations(iid_chr2, iid_lists)

# Optional: print summary
sapply(combo_table, length)

# View one combo
combo_table$type1_type2 

write.table(
  stack(combo_table),
  file = "genotype_combinations.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE,
  col.names = c("individual", "genotype_combo")
)

phenotypes <- read.delim("phenotype.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)


# First, assign type labels for Chr2
get_chr_type <- function(iid, iid_list) {
  for (i in 1:3) {
    if (iid %in% iid_list[[i]]) return(paste0("chr2_type", i))
  }
  return(NA)
}

# Repeat for Chr3
get_chr_type3 <- function(iid, iid_list) {
  for (i in 1:3) {
    if (iid %in% iid_list[[i]]) return(paste0("chr3_type", i))
  }
  return(NA)
}

# Apply the functions
phenotypes$chr2_type <- sapply(phenotypes$IID, get_chr_type, iid_list = iid_chr2)
phenotypes$chr3_type <- sapply(phenotypes$IID, get_chr_type3, iid_list = iid_chr3)

# Combine into a 9-genotype category
phenotypes$genotype_combo <- paste(phenotypes$chr2_type, phenotypes$chr3_type, sep = "_")
# Keep only rows where both chromosome types are assigned
phenotypes_clean <- subset(phenotypes, !is.na(chr2_type) & !is.na(chr3_type))
write.table(phenotypes_clean, file = "phenotypes_clean.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)

ggplot(phenotypes_clean, aes(x = genotype_combo, y = relMRWS)) +
  geom_boxplot(fill = "lightblue") +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(x = "Genotype Combination", y = "Relative MRWS") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(phenotypes_clean, aes(x = genotype_combo, y = relMRWS)) +
  geom_violin(fill = "lightgreen") +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(x = "Genotype Combination", y = "Relative MRWS") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(phenotypes_clean, aes(x = genotype_combo, y = relMRWS)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.6) +
  labs(x = "Genotype Combination", y = "Relative MRWS") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(phenotypes_clean, aes(x = chr3_type, y = relMRWS)) +
  geom_boxplot(fill = "lightblue") +
  facet_wrap(~ chr2_type) +
  labs(x = "Chr3 Genotype", y = "Relative MRWS", title = "Effect of Inversion Genotype on Chr2 and Chr3") +
  theme_minimal()

ggplot(phenotypes_clean, aes(x = chr3_type, y = relMRWS)) +
  geom_boxplot(fill = "lightblue") +
  facet_wrap(~ chr2_type) +
  labs(x = "Chr3 Genotype", y = "Relative MRWS", title = "Effect of Inversion Genotype on Chr2 and Chr3") +
  theme_minimal()

library(dplyr)

summary_stats <- phenotypes_clean %>%
  group_by(chr2_type, chr3_type) %>%
  summarise(
    mean_relMRWS = mean(relMRWS, na.rm = TRUE),
    sd = sd(relMRWS, na.rm = TRUE),
    n = n(),
    se = sd / sqrt(n)
  )

ggplot(summary_stats, aes(x = chr3_type, y = mean_relMRWS, group = chr2_type, color = chr2_type)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_relMRWS - se, ymax = mean_relMRWS + se), width = 0.2) +
  labs(x = "Chr3 Genotype", y = "Mean Relative MRWS", color = "Chr2 Genotype") +
  theme_minimal()

library(plotly)

plot_ly(
  data = plot_data,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~karyotype,
  colors = "Set1",  # or another palette
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4)
) %>%
  layout(
    title = "3D PCA Plot by Karyotype",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )



ggplot(plot_data, aes(x = PC1, y = PC2)) +
  geom_point(color = "black", size = 2) + 
  stat_ellipse(aes(color = karyotype), type = "norm", level = 0.95, linetype = "dashed", size = 1) +  # colored ellipses
  labs(
    title = "PCA with Colored Ellipses by Karyotype",
    x = "Principal Component 1 (PC1)",
    y = "Principal Component 2 (PC2)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
