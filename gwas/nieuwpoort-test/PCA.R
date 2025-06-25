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
