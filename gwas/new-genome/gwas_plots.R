####
# manhattan plots from lmm (using p_wald) or bslmm (bayesian, using PIP)

library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gridExtra)   # Arranging multiple grid-based plots
library(cowplot) # can align multi plots with plot_grid()
library(SNPRelate)
library(car)         # ANOVA and mixed models for factorial experiments
library(nlme)        # gls
library(emmeans)     # Estimated marginal means (least-squares means)
library(effects)     # Effect plots for regression models
library(scales)      # for label_percent()
library(htmlwidgets)
library(plotly)
library(car)        # ANOVA and mixed models for factorial experiments
library(tidyr)

#calculate chromosome coordinates
chrom.coords <- function(scafL,chromNames,gap = 9000000) {
  chromosome = vector()
  chromLengths  = vector()
  chromStarts = vector()
  chromEnds = vector()
  chromMid = vector()
  chrom = 1
  endLast = 0
  scafCurrent <- subset(scafL, chromosome == chromNames[1])
  chromosome[chrom] <- chrom
  chromLengths[chrom] <- sum(scafCurrent$length)
  chromStarts[chrom] <- endLast + 1
  chromEnds[chrom] <- endLast + chromLengths[chrom]
  chromMid[chrom] <- endLast + chromLengths[chrom]/2
  endLast = chromEnds[chrom]
  chrom = chrom + 1
  for (i in 2:length(chromNames)) {
    chromosome[chrom] <- chrom
    scafCurrent <- subset(scafL, chromosome == chromNames[i])
    chromLengths[chrom] <- sum(scafCurrent$length)
    chromStarts[chrom] <- endLast + gap + 1
    chromEnds[chrom] <- endLast + gap + chromLengths[chrom]
    chromMid[chrom] <- endLast + gap + chromLengths[chrom]/2
    endLast = chromEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(chromosome,chromLengths,chromStarts,chromEnds,chromMid))
  return(table)
}

#calculate scaffold coordinates
scaf.coords <- function(scafL,gap = 0) {
  scaffold = vector()
  scafStarts = vector()
  scafEnds = vector()
  chrom = 1
  endLast = 0
  for (e in 1:nrow(scafL)) {
    scaffold[chrom] <- levels(scafL$scaffold)[e]
    scafStarts[chrom] <- endLast + gap + 1
    scafEnds[chrom] <- endLast + gap + scafL$length[e]
    endLast = scafEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(scaffold,scafStarts,scafEnds))
  return(table)
}

# Load scaffold lengths
setwd("/scratch/leuven/357/vsc35707/popgen/new-genome-july-2025/")
scafL <- read.table("Pchalceus_SW.sorted.lengths.txt", header = FALSE, stringsAsFactors = FALSE)

# Add column names
colnames(scafL) <- c("chrom_name", "length")

# Extract CHR entries only
scafL <- scafL[grepl("^CHR", scafL$chrom_name), ]

# Make numeric
scafL$chrom_id <- as.numeric(gsub("CHR", "", scafL$chrom_name))

# Prepare input for function
scafL$chromosome <- scafL$chrom_id

# Make numeric
scafL$length <- as.numeric(scafL$length)

# Run function
chrom_coords <- chrom.coords(scafL = scafL, chromNames = 1:11)
print(chrom_coords)

# Merge chromosome and scaffold datasets
scaf_coords <- merge(scafL, chrom_coords, by="chromosome", all.x=TRUE)
names(scaf_coords)[names(scaf_coords) == 'chrom_name'] <- 'scaffold'
scaf_coords$scaffold<-as.character(scaf_coords$scaffold)

scaf_coords$chromosome <- as.factor(scaf_coords$chromosome)
#write.table(scaf_coords, file = "scaf_coords.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


## First plot the Fst scans for Spanish populations
setwd("/scratch/leuven/357/vsc35707/popgen/new-genome-july-2025/final-stats")
pattern_SP <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop3\\.stats"
file_list_SP <- list.files(pattern = pattern_SP, full.names = TRUE)
stats_SP <- rbindlist(lapply(file_list_SP, fread))
stats_SP <- as.data.frame(stats_SP)

comp <- list(stats_SP)
names <- c("Spain T/S")
col=c("black")

begin = chrom_coords$chromStarts[1]/1000000
end = max(chrom_coords$chromEnds)/1000000

# Set top and bottom plot limits
top <- 1
bot <- 0

# Number of chromosomes
num_chr <- 11

# Compute chromosome midpoints for labeling
chrom_coords$chromMids <- (chrom_coords$chromStarts + chrom_coords$chromEnds) / 2

# Color function for points
colfuncR <- colorRampPalette(c("black", "red"))
colL <- colfuncR(100)

# Set outer margins
par(mai = c(0.05, 0.4, 0.05, 0.4), oma = c(2, 0, 1, 0))

# Loop over all comparisons
for (i in 1:length(comp)) {
  # Merge comparison data with chromosome coordinates
  comb <- merge(comp[[i]], scaf_coords, by = 'scaffold', all.x = TRUE)
  comb <- as.data.frame(comb)
  comb$chromPos <- comb$mid + comb$chromStarts
  comb <- na.omit(comb)
  
  # Clean up the Fst column (assumed to be the 16th)
  comb <- comb %>% 
    filter(comb[, 16] != "none") %>%
    mutate(across(16, as.numeric))
  
  comb[, 16][is.na(comb[, 16])] <- 0
  comb[, 16][comb[, 16] < 0] <- 0
  comb$col <- colfuncR(100)[as.integer((comb[, 16] / 1) * 100) + 1]
  
  # Set up an empty plot
  plot(0, pch = "", xlim = c(min(chrom_coords$chromStarts), max(chrom_coords$chromEnds)) / 1e6,
       ylim = c(bot, top), ylab = "", yaxt = "n", lwd = 0.5,
       xlab = "", xaxt = "n", bty = "n", main = "", axes = FALSE)
  
  # Add alternating gray background per chromosome
  rect(chrom_coords$chromStarts / 1e6, rep(bot, num_chr),
       chrom_coords$chromEnds / 1e6, rep(top, num_chr),
       col = rep(c("gray95", "gray90"), length.out = num_chr),
       border = NA)
  
  # Overlay Fst points
  par(new = TRUE)
  plot(comb$chromPos / 1e6, comb[, 16], type = "p", pch = 19, cex = 0.7,
       col = adjustcolor(comb$col, alpha.f = 1),
       xlim = c(min(chrom_coords$chromStarts), max(chrom_coords$chromEnds)) / 1e6,
       ylim = c(bot, top), axes = FALSE, bty = "n", xlab = "", ylab = "")
  
  # Add y-axis
  axis(2, cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(F[ST]), cex = 0.5, line = 0.8)
  
  # Add comparison label
  mtext(side = 4, text = names[i], cex = 0.6)
}

# Add chromosome number labels on x-axis
axis(1, at = chrom_coords$chromMids[1:num_chr] / 1e6,
     labels = 1:num_chr, lwd = 0, lwd.ticks = 0)

# Add scale bar
segments(x0 = 5, y0 = 0.93, x1 = 15, y1 = 0.93, lwd = 1)
text(10, 0.95, labels = "10Mb", cex = 1)


comb
#write.table(comb, file = "spain_stats.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

chrom_colors <- rep(c("lightgray", "darkgray"), length.out = n_chrom_1)

Fst <- ggplot(comb,aes(x = chromPos,y = FstWC_S_HUE_T_S_COT_S_egg)) +
  geom_point(aes(color = factor(chromosome))) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  scale_y_continuous(
    breaks = function(x) {
      brks <- pretty(x)
      brks[seq(1, length(brks), by = 2)]  
    }
  )+
  theme_minimal() +
  labs(x = "Chromosome",
       y = expression(F[st])) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16))
Fst

# The make the Manhattan plots
gemma_results_NP25_1 <- fread("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/gemma_results/gemma_lmm_results.assoc.txt")

gemma_results_NP25_1$chr <- as.factor(gemma_results_NP25_1$chr)  # Make sure chromosome is a factor
gemma_results_NP25_1$p_wald <- as.numeric(gemma_results_NP25_1$p_wald)
gemma_results_NP25_1$ps <- as.numeric(gemma_results_NP25_1$ps)

gemma_results_NP25_1 <- gemma_results_NP25_1 %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)

write.table(gemma_results_NP25_1, file = "gemma_results_NP25_MRWS.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# bonferroni treshold corrects for number of tests (conservative)
num_tests_1 <- nrow(gemma_results_NP25_1)  # Total number of SNPs tested
bonferroni_threshold_1 <- 0.05 / num_tests_1  # Genome-wide significance threshold
bonferroni_threshold_1_log <- -log10(bonferroni_threshold_1)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_1 <- length(unique(gemma_results_NP25_1$chr))
chrom_colors <- rep(c("lightgray", "darkgray"), length.out = n_chrom_1)

# Plot1
manhattan_pwald_1 <- ggplot(gemma_results_NP25_1, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_1, p_wald > bonferroni_threshold_1),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_1, p_wald < bonferroni_threshold_1),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_1_log, color = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", size = 0.5) +  
  annotate("text", x = max(gemma_results_NP25_1$pos_cumulative) * 0.88, 
           y = bonferroni_threshold_1_log + .5, label = "Bonferroni threshold", color = "red", size = 4) +
  annotate("text", x = max(gemma_results_NP25_1$pos_cumulative) * 0.873, 
           y = -log10(1e-5) - .5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  annotate("segment", x = 5e6, xend = 25e6, y = 18, yend = 18, size = 0.5) +
  annotate("text", x = 15e6, y = 17, label = "20 Mb", size = 3.5) +
  theme_minimal() +
  labs(title = "%MRWS (n=192)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_pwald_1



manhattanANDfst <- plot_grid(
  manhattan_pwald_1,
  Fst,
  ncol = 1,
  align = 'v',
  axis = 'lr',
  rel_heights = c(1, 0.4)
)
manhattanANDfst

gwas_wl = subset(gemma_results_NP25_1, p_wald < bonferroni_threshold_1*100)
gwas_wl_chr5 = subset(gwas_wl, chr == "5")
pos_min = min(gwas_wl_chr5$ps)
pos_max = max(gwas_wl_chr5$ps)
pos_min
pos_max

gwas_wl_chr6 = subset(gwas_wl, chr == "6")
pos_min = min(gwas_wl_chr6$ps)
pos_max = max(gwas_wl_chr6$ps)
pos_min
pos_max

setwd("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/gemma_results/")
# Load param.txt, %MRWS n=192
param <- fread("gemma_bslmm_results.param.txt")

colnames(param) <- c("chr", "SNP", "ps", "n_miss", "alpha", "beta", "PIP")
param$chr <- as.factor(param$chr)

# Optional: add cumulative position for Manhattan plot
param <- param %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)

# Plot BSLMM, %MRWS n=192
manhattan_PIP <- ggplot(param, aes(x = pos_cumulative, y = PIP)) +
  geom_point(data = subset(param, PIP < 0.1),
             aes(x = pos_cumulative, y = PIP, color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(param, PIP > 0.1),
             aes(x = pos_cumulative, y = PIP),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed", size = 0.5) +  
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  scale_y_sqrt() + # square root transformation of the y axis
  theme_minimal() +
  labs(title = "BSLMM: %MRWS (n=192)",
       x = NULL,
       y = expression(sqrt(PIP))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_PIP

# %MRWS chr5:AA n=120
param <- fread("gemma_bslmm_results_chr5_homoAA.param.txt")
colnames(param) <- c("chr", "SNP", "ps", "n_miss", "alpha", "beta", "PIP")
param$chr <- as.factor(param$chr)
param <- param %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)

# Plot BSLMM, %MRWS chr5:AA n=120
manhattan_PIP_AA <- ggplot(param, aes(x = pos_cumulative, y = PIP)) +
  geom_point(data = subset(param, PIP < 0.1),
             aes(x = pos_cumulative, y = PIP, color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(param, PIP > 0.1),
             aes(x = pos_cumulative, y = PIP),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed", size = 0.5) +  
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  scale_y_sqrt() + # square root transformation of the y axis
  theme_minimal() +
  labs(title = "BSLMM: %MRWS - chr5:AA (n=120)",
       x = NULL,
       y = expression(sqrt(PIP))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_PIP_AA

# EL, n=192
param <- fread("gemma_bslmm_results_EL.param.txt")
colnames(param) <- c("chr", "SNP", "ps", "n_miss", "alpha", "beta", "PIP")
param$chr <- as.factor(param$chr)
param <- param %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)
manhattan_PIP_EL <- ggplot(param, aes(x = pos_cumulative, y = PIP)) +
  geom_point(data = subset(param, PIP < 0.1),
             aes(x = pos_cumulative, y = PIP, color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(param, PIP > 0.1),
             aes(x = pos_cumulative, y = PIP),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed", size = 0.5) +  
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  scale_y_sqrt() + # square root transformation of the y axis
  theme_minimal() +
  labs(title = "BSLMM: Elytron length (n=192)",
       x = NULL,
       y = expression(sqrt(PIP))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_PIP_EL


# Average emergence time, n=192
param <- fread("gemma_bslmm_results_time.param.txt")
colnames(param) <- c("chr", "SNP", "ps", "n_miss", "alpha", "beta", "PIP")
param$chr <- as.factor(param$chr)
param <- param %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)
manhattan_PIP_time <- ggplot(param, aes(x = pos_cumulative, y = PIP)) +
  geom_point(data = subset(param, PIP < 0.1),
             aes(x = pos_cumulative, y = PIP, color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(param, PIP > 0.1),
             aes(x = pos_cumulative, y = PIP),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = 0.1, color = "red", linetype = "dashed", size = 0.5) +  
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  scale_y_sqrt() + # square root transformation of the y axis
  theme_minimal() +
  labs(title = "BSLMM: Average emergence time (n=192, 3 trials)",
       x = NULL,
       y = expression(sqrt(PIP))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_PIP_time


setwd("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/pos/")

# chr5
vcfFile_chr5 <- "Pchal_Bar_SW.chr_5.window50Mb.vcf.gz"
gdsFile_chr5 <- "Pchal_Bar_SW.chr_5.window50Mb.gds"

# Convert VCF to GDS format
snpgdsVCF2GDS(vcfFile_chr5, gdsFile_chr5, method = "biallelic.only")

# Open GDS file for analysis
genofile_chr5 <- snpgdsOpen(gdsFile_chr5)

# Perform PCA analysis on the genotype data
ccm_pca <- snpgdsPCA(genofile_chr5, autosome.only = TRUE)

# Extract sample names and clean them by removing extra suffixes (whatever you have in your sample names, you maybe cleaned them already)
sNames <- ccm_pca$sample.id
sNames <- sub(".Pchal_Bar_SW.filtered.sorted.bam", "", sNames)
sNames <- sub("bams/", "", sNames)

popInfo <- data.frame(
  V1 = c("NP25")
)

# Create a data frame for plotting PCA results
plot_data <- as.data.frame(ccm_pca$eigenvect)
colnames(plot_data) <- paste0("PC", 1:ncol(plot_data))  # Rename columns for clarity

# divide in three clusters based on PC1
plot_data$karyotype <- cut(plot_data$PC1, breaks = 3, labels = c("AA", "Aa", "aa"))
plot_data$IID <- ccm_pca$sample.id

# PCA Plot Colored by Population
ggplot(plot_data, aes(x = PC1, y = PC2, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Karyotype, 20 Mb window in Chromosome 5", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# PC2 and PC3
ggplot(plot_data, aes(x = PC2, y = PC3, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Karyotype, 20 Mb window in Chromosome 5", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# IID lists of different karyotypes
iid_lists <- split(plot_data$IID, plot_data$karyotype)
iid_lists$AA
iid_lists$Aa
iid_lists$aa

#write.table(iid_lists$AA, "homozygotes_chr5_AA.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(iid_lists$Aa, "heterozygotes_chr5_Aa.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(iid_lists$aa, "homozygotes_chr5_aa.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# chr6
vcfFile_chr6 <- "Pchal_Bar_SW.chr_6.window15kb.vcf.gz"
gdsFile_chr6_15Kb <- "Pchal_Bar_SW.chr_6.window15kb.gds"

# Convert VCF to GDS format
snpgdsVCF2GDS(vcfFile_chr6, gdsFile_chr6_15Kb, method = "biallelic.only")

# Open GDS file for analysis
genofile_chr6 <- snpgdsOpen(gdsFile_chr6_15Kb)

# Perform PCA analysis on the genotype data
ccm_pca <- snpgdsPCA(genofile_chr6, autosome.only = TRUE)

# Extract sample names and clean them by removing extra suffixes (whatever you have in your sample names, you maybe cleaned them already)
sNames <- ccm_pca$sample.id
sNames <- sub(".Pchal_Bar_SW.filtered.sorted.bam", "", sNames)
sNames <- sub("bams/", "", sNames)

popInfo <- data.frame(
  V1 = c("NP25")
)

# Create a data frame for plotting PCA results
plot_data <- as.data.frame(ccm_pca$eigenvect)
colnames(plot_data) <- paste0("PC", 1:ncol(plot_data))  # Rename columns for clarity

# divide in three clusters based on PC1
plot_data$karyotype <- cut(plot_data$PC1, breaks = 3, labels = c("AA", "Aa", "aa"))
plot_data$IID <- ccm_pca$sample.id

# PCA Plot Colored by Population
ggplot(plot_data, aes(x = PC1, y = PC2, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Karyotype, 15 Kb window in Chromosome 6", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# PC2 and PC3
ggplot(plot_data, aes(x = PC2, y = PC3, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Karyotype, 15 Kb window in Chromosome 6", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# IID lists of different karyotypes
iid_lists <- split(plot_data$IID, plot_data$karyotype)
iid_lists$AA
iid_lists$Aa
iid_lists$aa

#write.table(iid_lists$AA, "homozygotes_chr6_AA.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(iid_lists$Aa, "heterozygotes_chr6_Aa.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(iid_lists$aa, "homozygotes_chr6_aa.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Observed genotype counts
obs <- c(AA = 11, Aa = 61, aa = 120)
n <- sum(obs)

# Calculate allele frequencies
total_alleles <- 2 * n
p <- (2 * obs["AA"] + obs["Aa"]) / total_alleles  # frequency of A
q <- 1 - p                                        # frequency of a

# Expected genotype frequencies under HWE
exp <- c(
  AA = p^2 * n,
  Aa = 2 * p * q * n,
  aa = q^2 * n
)

# Chi-square test
chisq_result <- chisq.test(x = obs, p = exp / n, rescale.p = TRUE)

# Output
cat("Observed counts:\n")
print(obs)
cat("\nExpected counts under HWE:\n")
print(round(exp, 2))
cat("\nChi-square test result:\n")
print(chisq_result)



# The make the Manhattan plots
gemma_results_NP25_2 <- fread("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/gemma_results/gemma_lmm_results_cov.assoc.txt")

gemma_results_NP25_2$chr <- as.factor(gemma_results_NP25_2$chr)  # Make sure chromosome is a factor
gemma_results_NP25_2$p_wald <- as.numeric(gemma_results_NP25_2$p_wald)
gemma_results_NP25_2$ps <- as.numeric(gemma_results_NP25_2$ps)

gemma_results_NP25_2 <- gemma_results_NP25_2 %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)


# bonferroni treshold corrects for number of tests (conservative)
num_tests_2 <- nrow(gemma_results_NP25_2)  # Total number of SNPs tested
bonferroni_threshold_2 <- 0.05 / num_tests_2  # Genome-wide significance threshold
bonferroni_threshold_2_log <- -log10(bonferroni_threshold_2)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_2 <- length(unique(gemma_results_NP25_2$chr))
n_chrom_2 <- rep(c("lightgray", "darkgray"), length.out = n_chrom_2)

# Plot1
manhattan_pwald_2 <- ggplot(gemma_results_NP25_2, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_2, p_wald > bonferroni_threshold_2),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_2, p_wald < bonferroni_threshold_2),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_2_log, color = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", size = 0.5) +  
  annotate("text", x = max(gemma_results_NP25_2$pos_cumulative) * 0.88, 
           y = bonferroni_threshold_2_log + .5, label = "Bonferroni threshold", color = "red", size = 4) +
  annotate("text", x = max(gemma_results_NP25_2$pos_cumulative) * 0.873, 
           y = -log10(1e-5) - .5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  annotate("segment", x = 5e6, xend = 25e6, y = 18, yend = 18, size = 0.5) +
  annotate("text", x = 15e6, y = 17, label = "20 Mb", size = 3.5) +
  theme_minimal() +
  labs(title = "%MRWS (n=192)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_pwald_2
manhattan_pwald_1



# The make the Manhattan plots
gemma_results_NP25_3 <- fread("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/gemma_results/gemma_lmm_results_chr5_homoAA.assoc.txt")

gemma_results_NP25_3$chr <- as.factor(gemma_results_NP25_3$chr)  # Make sure chromosome is a factor
gemma_results_NP25_3$p_wald <- as.numeric(gemma_results_NP25_3$p_wald)
gemma_results_NP25_3$ps <- as.numeric(gemma_results_NP25_3$ps)

gemma_results_NP25_3 <- gemma_results_NP25_3 %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)

#write.table(gemma_results_NP25_3, file = "gemma_results_NP25_chr5AA.tsv", sep = "\t", quote = FALSE, row.names = FALSE)



# bonferroni treshold corrects for number of tests (conservative)
num_tests_3 <- nrow(gemma_results_NP25_3)  # Total number of SNPs tested
bonferroni_threshold_3 <- 0.05 / num_tests_3  # Genome-wide significance threshold
bonferroni_threshold_3_log <- -log10(bonferroni_threshold_3)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_3 <- length(unique(gemma_results_NP25_3$chr))
n_chrom_3 <- rep(c("lightgray", "darkgray"), length.out = n_chrom_3)

# Plot1
manhattan_pwald_3 <- ggplot(gemma_results_NP25_3, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_3, p_wald > bonferroni_threshold_3),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_3, p_wald < bonferroni_threshold_3),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_3_log, color = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", size = 0.5) +  
  annotate("text", x = max(gemma_results_NP25_3$pos_cumulative) * 0.88, 
           y = bonferroni_threshold_3_log + .5, label = "Bonferroni threshold", color = "red", size = 4) +
  annotate("text", x = max(gemma_results_NP25_3$pos_cumulative) * 0.873, 
           y = -log10(1e-5) - .5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  annotate("segment", x = 5e6, xend = 25e6, y = 18, yend = 18, size = 0.5) +
  annotate("text", x = 15e6, y = 17, label = "20 Mb", size = 3.5) +
  theme_minimal() +
  labs(title = "%MRWS - chr5:AA individuals (n=120)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_pwald_3
manhattan_pwald_2
manhattan_pwald_1


manhattanANDfst_2 <- plot_grid(
  manhattan_pwald_2,
  manhattan_pwald_3,
  Fst,
  ncol = 1,
  align = 'v',
  axis = 'lr',
  rel_heights = c(0.8, 0.8, 0.6)
)
manhattanANDfst_2


gwas_wl_3 = subset(gemma_results_NP25_3, p_wald < bonferroni_threshold_3*100)
gwas_wl_chr4 = subset(gwas_wl_3, chr == "4")
pos_min = min(gwas_wl_chr4$ps)
pos_max = max(gwas_wl_chr4$ps)
pos_min
pos_max

setwd("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/pos/")
# chr4
# use lead significant SNPs (on positions 33865711 and 33872519) to genotype these individuals
vcfFile_chr4_snp2 <- "Pchal_Bar_SW.chr_4.snp2.vcf.gz"
gdsFile_chr4_snp2 <- "Pchal_Bar_SW.chr_4.snp2.gds"

# Convert VCF to GDS format
snpgdsVCF2GDS(vcfFile_chr4_snp2, gdsFile_chr4_snp2, method = "biallelic.only")

# Open GDS file for analysis
genofile_chr4 <- snpgdsOpen(gdsFile_chr4_snp2)

# Perform PCA analysis on the genotype data
ccm_pca <- snpgdsPCA(genofile_chr4, autosome.only = TRUE)

# Extract sample names and clean them by removing extra suffixes (whatever you have in your sample names, you maybe cleaned them already)
sNames <- ccm_pca$sample.id
sNames <- sub(".Pchal_Bar_SW.filtered.sorted.bam", "", sNames)
sNames <- sub("bams/", "", sNames)

popInfo <- data.frame(
  V1 = c("NP25")
)

# Create a data frame for plotting PCA results
plot_data <- as.data.frame(ccm_pca$eigenvect)
colnames(plot_data) <- paste0("PC", 1:ncol(plot_data))  # Rename columns for clarity

# divide in three clusters based on PC1
plot_data$karyotype <- cut(plot_data$PC1, breaks = 3, labels = c("BB", "Bb", "bb"))
plot_data$IID <- ccm_pca$sample.id

# PCA Plot Colored by Population
ggplot(plot_data, aes(x = PC1, y = PC2, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Karyotype, SNP 33865711 in Chromosome 4", x = "Principal Component 1 (PC1)", y = "Principal Component 2 (PC2)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# PC2 and PC3
ggplot(plot_data, aes(x = PC2, y = PC3, color = karyotype)) +
  geom_point(size = 2) +
  labs(title = "PCA Colored by Karyotype, SNP 33865711 in Chromosome 4", x = "Principal Component 2 (PC2)", y = "Principal Component 3 (PC3)") +
  theme_minimal() +
  theme(legend.title = element_blank())

# IID lists of different karyotypes
iid_lists <- split(plot_data$IID, plot_data$karyotype)
iid_lists$BB
iid_lists$Bb
iid_lists$bb

#write.table(iid_lists$BB, "homozygotes_chr4_BB.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(iid_lists$Bb, "heterozygotes_chr4_Bb.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)
#write.table(iid_lists$bb, "homozygotes_chr4_bb.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)


fig <- plot_ly(
  data = plot_data,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~karyotype,
  colors = c("red", "green", "blue"),  # AA, Aa, aa (adjust as needed)
  text = ~IID,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4)
)

# Add axis labels
fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "PC1"),
    yaxis = list(title = "PC2"),
    zaxis = list(title = "PC3")
  )
)

fig

# Observed genotype counts
obs <- c(BB = 12, Bb = 85, bb = 95)
n <- sum(obs)

# Calculate allele frequencies
total_alleles <- 2 * n
p <- (2 * obs["BB"] + obs["Bb"]) / total_alleles  # frequency of B
q <- 1 - p                                        # frequency of b

# Expected genotype frequencies under HWE
exp <- c(
  BB = p^2 * n,
  Bb = 2 * p * q * n,
  bb = q^2 * n
)

# Chi-square test
chisq_result <- chisq.test(x = obs, p = exp / n, rescale.p = TRUE)

# Output
cat("Observed counts:\n")
print(obs)
cat("\nExpected counts under HWE:\n")
print(round(exp, 2))
cat("\nChi-square test result:\n")
print(chisq_result)


## Test for epistasis between chr5 and chr4
setwd("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort")
phenotypes <- fread("phenotype_renamed.txt")
phenotypes <- na.omit(phenotypes)

setwd("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/pos")
# chr5
genotype_chr5_aa1 <- fread("homozygotes_chr5_aa.txt", header = FALSE)
genotype_chr5_aA2 <- fread("heterozygotes_chr5_Aa.txt", header = FALSE)
genotype_chr5_AA3 <- fread("homozygotes_chr5_AA.txt", header = FALSE)

#chr4
genotype_chr4_bb1 <- fread("homozygotes_chr4_bb.txt", header = FALSE)
genotype_chr4_bB2 <- fread("heterozygotes_chr4_Bb.txt", header = FALSE)
genotype_chr4_BB3 <- fread("homozygotes_chr4_BB.txt", header = FALSE)


phenotypes$chr5 <- NA
phenotypes$chr5[phenotypes$IID %in% genotype_chr5_aa1$V1] <- "chr5:aa"
phenotypes$chr5[phenotypes$IID %in% genotype_chr5_aA2$V1] <- "chr5:aA"
phenotypes$chr5[phenotypes$IID %in% genotype_chr5_AA3$V1] <- "chr5:AA"

phenotypes$chr4 <- NA
phenotypes$chr4[phenotypes$IID %in% genotype_chr4_bb1$V1] <- "chr4:bb"
phenotypes$chr4[phenotypes$IID %in% genotype_chr4_bB2$V1] <- "chr4:bB"
phenotypes$chr4[phenotypes$IID %in% genotype_chr4_BB3$V1] <- "chr4:BB"

# Combine chr5 and chr4 into a single factor
phenotypes$genotype <- interaction(phenotypes$chr5, phenotypes$chr4, drop = TRUE)

# Clean up labels (e.g., "aA.AA" → "aA x AA")
levels(phenotypes$genotype) <- gsub("\\.", "_x_", levels(phenotypes$genotype))

phenotypes$genotype <- factor(
  phenotypes$genotype,
  levels = c("chr5:aa_x_chr4:bb", "chr5:aa_x_chr4:bB", "chr5:aa_x_chr4:BB",
             "chr5:aA_x_chr4:bb", "chr5:aA_x_chr4:bB", "chr5:aA_x_chr4:BB",
             "chr5:AA_x_chr4:bb", "chr5:AA_x_chr4:bB", "chr5:AA_x_chr4:BB"))

plot_epistasis <- ggplot(phenotypes, aes(x = genotype, y = relMRWS)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  labs(
    x = "Genotype (chr5 x chr4)",
    y = "%MRWS",
    title = "%MRWS across genotype combinations"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major.x = element_blank()
  )
plot_epistasis

fit <- lm(relMRWS ~ genotype, data = phenotypes)
Anova(fit, type='III')
ncvTest(fit) # non equal variance!

fit <- gls(relMRWS ~ genotype, data = phenotypes)
Anova(fit, type='III')


# Make sure genotype columns are factors without the "chrX:" prefix
phenotypes$chr5_geno <- factor(sub("chr5:", "", phenotypes$chr5),
                               levels = c("aa", "aA", "AA"))
phenotypes$chr4_geno <- factor(sub("chr4:", "", phenotypes$chr4),
                               levels = c("bb", "bB", "BB"))

# Check the structure
table(phenotypes$chr5_geno, phenotypes$chr4_geno)

# Fit GLS model with interaction
fit_gls <- gls(relMRWS ~ chr5_geno * chr4_geno,
               data = phenotypes)

# Type III ANOVA (interaction test)
Anova(fit_gls, type = "III")
summary(fit_gls)

# Model with heteroscedasticity structure. This allows different variances per genotype combination
# Make a factor for the combination of genotypes
phenotypes$geno_group <- interaction(phenotypes$chr5_geno, phenotypes$chr4_geno, drop = TRUE)

# Fit GLS model allowing different variances for each genotype combination
fit_gls_var <- gls(
  relMRWS ~ chr5_geno * chr4_geno,
  data = phenotypes,
  weights = varIdent(form = ~ 1 | geno_group)
)

# Type III ANOVA
Anova(fit_gls_var, type = "III")

# Compare AICs
AIC(fit_gls, fit_gls_var)

# Compute means and counts together
counts_means_df <- phenotypes %>%
  group_by(chr5_geno, chr4_geno) %>%
  summarise(mean_relMRWS = mean(relMRWS, na.rm = TRUE),
            n = n(), .groups = "drop") %>%
  mutate(label = paste0("n = ", n))

ggplot(phenotypes, aes(x = chr5_geno, y = relMRWS,
                       color = chr4_geno, group = chr4_geno)) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(width = 0.3)) +
  stat_summary(fun = mean, geom = "line", position = position_dodge(width = 0.3)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2,
               position = position_dodge(width = 0.3)) +
  # Add counts above each mean
  #geom_text(data = counts_means_df,
  #          aes(x = chr5_geno,
  #              y = mean_relMRWS + 0.0025,  # slightly above the mean
  #              label = label,
  #              color = chr4_geno),
  #          position = position_dodge(width = 0.2),
  #          inherit.aes = FALSE) +
  theme_minimal() +
  labs(x = "Genotype at chr5 SNP",
       y = "Relative MRWS",
       color = "Genotype at chr4 SNP") +
  ggtitle("Interaction between chr5 and chr4 SNPs")

# formal test of epistasis (Sander's code), using GLS instead of LM 
fit1 <- gls(relMRWS ~ chr5 * chr4, data=phenotypes)
AIC(fit1)
fit2 <- gls(relMRWS ~ chr5 + chr4, data=phenotypes) # best model has the additive effect, not interaction
AIC(fit2)
fit3 <- gls(relMRWS ~ chr4, data=phenotypes) # worst model, the SNP on chr4 has no effect on wing size on its own
AIC(fit3)
fit4 <- gls(relMRWS ~ chr5, data=phenotypes) # second best model, just chr5 
AIC(fit4)

Anova(fit1, type='III')
summary(fit1)
plot(allEffects(fit1))

Anova(fit2, type='III')
summary(fit2)
plot(allEffects(fit2))

Anova(fit3, type='III')
summary(fit3)
plot(allEffects(fit3))

Anova(fit4, type='III')
summary(fit4)
plot(allEffects(fit4))


table(phenotypes$chr4_geno, phenotypes$chr5_geno)
chisq.test(phenotypes$chr4_geno, phenotypes$chr5_geno)
fisher.test(table(phenotypes$chr4_geno, phenotypes$chr5_geno))

phenotypes$fitted <- fitted(fit1)

ggplot(phenotypes, aes(x = chr5, y = fitted, color = chr4, group = chr4)) +
  stat_summary(fun = mean, geom = "line", size = 1.2) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  labs(title = "Fitted Interaction: chr5 × chr4", y = "Fitted relMRWS") +
  theme_minimal()


emm <- emmeans(fit1, ~ chr4 * chr5)
emm
counts <- phenotypes %>%
  group_by(chr4, chr5) %>%
  summarise(n = n(), .groups = "drop")
counts <- as.data.frame(counts)
emm <- as.data.frame(emm)

# Convert emm to a data frame
emm_df <- as.data.frame(emm)

# Merge on chr4 and chr5
emm_df <- merge(emm, counts, by = c("chr4", "chr5"), all.x = TRUE)

emm_df

# 5. Plot with sample sizes as labels
relMRWSchr5chr4 <- ggplot(emm_df, aes(x = chr5, y = emmean, color = chr4, group = chr4)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.5), width = 0.2) +
  geom_line(position = position_dodge(width = 0.5)) +
  geom_text(aes(label = paste0("n=", n)),
            position = position_dodge(width = 0.5),
            vjust = -1.2, size = 3.5, show.legend = FALSE) +
  labs(
    y = "Fitted %MRWS",
    x = NULL,
    title = NULL) +
  theme_minimal(base_size = 18) +
  theme(legend.position = c(0.6,0.9),legend.title = element_blank())
relMRWSchr5chr4



## GWAS for EL
# The make the Manhattan plots
gemma_results_NP25_EL <- fread("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/gemma_results/gemma_lmm_results_EL.assoc.txt")

gemma_results_NP25_EL$chr <- as.factor(gemma_results_NP25_EL$chr)  # Make sure chromosome is a factor
gemma_results_NP25_EL$p_wald <- as.numeric(gemma_results_NP25_EL$p_wald)
gemma_results_NP25_EL$ps <- as.numeric(gemma_results_NP25_EL$ps)

gemma_results_NP25_EL <- gemma_results_NP25_EL %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)


# bonferroni treshold corrects for number of tests (conservative)
num_tests_EL <- nrow(gemma_results_NP25_EL)  # Total number of SNPs tested
bonferroni_threshold_EL <- 0.05 / num_tests_EL  # Genome-wide significance threshold
bonferroni_threshold_EL_log <- -log10(bonferroni_threshold_EL)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_1 <- length(unique(gemma_results_NP25_EL$chr))
chrom_colors <- rep(c("lightgray", "darkgray"), length.out = n_chrom_1)

# Plot1
manhattan_pwald_EL <- ggplot(gemma_results_NP25_EL, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_EL, p_wald > bonferroni_threshold_EL),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_EL, p_wald < bonferroni_threshold_EL),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_EL_log, color = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", size = 0.5) +  
  annotate("text", x = max(gemma_results_NP25_EL$pos_cumulative) * 0.88, 
           y = bonferroni_threshold_EL_log + .5, label = "Bonferroni threshold", color = "red", size = 4) +
  annotate("text", x = max(gemma_results_NP25_EL$pos_cumulative) * 0.873, 
           y = -log10(1e-5) - .5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  #annotate("segment", x = 5e6, xend = 25e6, y = 21, yend = 21, size = 0.5) +
  #annotate("text", x = 15e6, y = 19.5, label = "20 Mb", size = 3.5) +
  theme_minimal() +
  labs(title = "Elytron length (n=192)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_pwald_EL

# Adding sex as a covariate
gemma_results_NP25_EL_sex <- fread("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/gemma_results/gemma_lmm_results_EL_sex.assoc.txt")
gemma_results_NP25_EL_sex$chr <- as.factor(gemma_results_NP25_EL_sex$chr)  # Make sure chromosome is a factor
gemma_results_NP25_EL_sex$p_wald <- as.numeric(gemma_results_NP25_EL_sex$p_wald)
gemma_results_NP25_EL_sex$ps <- as.numeric(gemma_results_NP25_EL_sex$ps)
gemma_results_NP25_EL_sex <- gemma_results_NP25_EL_sex %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)


# bonferroni treshold corrects for number of tests (conservative)
num_tests_EL_sex <- nrow(gemma_results_NP25_EL_sex)  # Total number of SNPs tested
bonferroni_threshold_EL_sex <- 0.05 / num_tests_EL_sex  # Genome-wide significance threshold
bonferroni_threshold_EL_log_sex <- -log10(bonferroni_threshold_EL_sex)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_1 <- length(unique(gemma_results_NP25_EL_sex$chr))
chrom_colors <- rep(c("lightgray", "darkgray"), length.out = n_chrom_1)

# Plot1
manhattan_pwald_EL_sex <- ggplot(gemma_results_NP25_EL_sex, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_EL_sex, p_wald > bonferroni_threshold_EL_sex),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_EL_sex, p_wald < bonferroni_threshold_EL_sex),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_EL_log_sex, color = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", size = 0.5) +  
  annotate("text", x = max(gemma_results_NP25_EL_sex$pos_cumulative) * 0.88, 
           y = bonferroni_threshold_EL_log_sex + .5, label = "Bonferroni threshold", color = "red", size = 4) +
  annotate("text", x = max(gemma_results_NP25_EL_sex$pos_cumulative) * 0.873, 
           y = -log10(1e-5) - .5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  #annotate("segment", x = 5e6, xend = 25e6, y = 7.5, yend = 7.5, size = 0.5) +
  #annotate("text", x = 15e6, y = 7, label = "20 Mb", size = 3.5) +
  theme_minimal() +
  labs(title = "Elytron length (n=192, sex as covariate)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_pwald_EL_sex


elytra <- plot_grid(
  manhattan_pwald_EL,
  manhattan_pwald_EL_sex,
  manhattan_PIP_EL,
  ncol = 1,
  align = 'v',
  axis = 'lr',
  rel_heights = c(1, 1, 1)
)
elytra

## GWAS for average emergence time
# The make the Manhattan plots
gemma_results_NP25_time <- fread("/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/gemma_results/gemma_lmm_results_time.assoc.txt")

gemma_results_NP25_time$chr <- as.factor(gemma_results_NP25_time$chr)  # Make sure chromosome is a factor
gemma_results_NP25_time$p_wald <- as.numeric(gemma_results_NP25_time$p_wald)
gemma_results_NP25_time$ps <- as.numeric(gemma_results_NP25_time$ps)

gemma_results_NP25_time <- gemma_results_NP25_time %>%
  left_join(scaf_coords, by = c("chr" = "chromosome")) %>%
  mutate(pos_cumulative = ps + chromStarts)


# bonferroni treshold corrects for number of tests (conservative)
num_tests_time <- nrow(gemma_results_NP25_time)  # Total number of SNPs tested
bonferroni_threshold_time  <- 0.05 / num_tests_time  # Genome-wide significance threshold
bonferroni_threshold_time_log <- -log10(bonferroni_threshold_time)  # Convert to -log10 scale

# alternating grayscale for adjacent chromosomes
n_chrom_1 <- length(unique(gemma_results_NP25_time$chr))
chrom_colors <- rep(c("lightgray", "darkgray"), length.out = n_chrom_1)

# Plot1
manhattan_pwald_time <- ggplot(gemma_results_NP25_time, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(data = subset(gemma_results_NP25_time, p_wald > bonferroni_threshold_time),
             aes(x = pos_cumulative, y = -log10(p_wald), color = factor(chr)),
             size = 2, alpha = 1) +
  geom_point(data = subset(gemma_results_NP25_time, p_wald < bonferroni_threshold_time),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 1) + # Red for significant SNPs
  geom_hline(yintercept = bonferroni_threshold_time_log, color = "red", linetype = "dashed", size = 0.5) + 
  geom_hline(yintercept = -log10(1e-5), color = "blue", linetype = "dashed", size = 0.5) +  
  annotate("text", x = max(gemma_results_NP25_time$pos_cumulative) * 0.88, 
           y = bonferroni_threshold_time_log + .5, label = "Bonferroni threshold", color = "red", size = 4) +
  annotate("text", x = max(gemma_results_NP25_time$pos_cumulative) * 0.873, 
           y = -log10(1e-5) - .5, label = "Suggestive threshold", color = "blue", size = 4) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = scaf_coords$chromMid,
    labels = scaf_coords$chromosome) +
  annotate("segment", x = 5e6, xend = 25e6, y = 7.5, yend = 7.5, size = 0.5) +
  annotate("text", x = 15e6, y = 7, label = "20 Mb", size = 3.5) +
  theme_minimal() +
  labs(title = "Average emergence time (n=192, 3 trials, sex as covariate)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_pwald_time
