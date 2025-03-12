#install.packages("data.table")

setwd("/vsc-hard-mounts/leuven-data/357/vsc35707")

library(data.table)

# Results for ALL data points
gemma_results <- fread("gemma_lmm_results.assoc.txt")

gemma_results$chr <- as.factor(gemma_results$chr)  # Make sure chromosome is a factor
gemma_results$p_wald <- as.numeric(gemma_results$p_wald)  # p-value

library(ggplot2)

# Set a threshold for significance (optional)
threshold <- 0.000001
threshold_2 <- 0.0001

# Create a data frame for chromosome labels and positions
chr_labels <- data.frame(chr = levels(gemma_results$chr),
                         label = paste("chr", levels(gemma_results$chr), sep = ""))

# Create the Manhattan plot
ggplot(gemma_results, aes(x = ps, y = -log10(p_wald))) +
  geom_point(aes(color = chr), size = 1, alpha = 0.6) +
  scale_color_manual(values = rep("blue", length(unique(gemma_results$chr)))) +  # Adjust color if needed
  theme_minimal() +
  labs(title = "Manhattan Plot",
       x = "Chromosome Position",
       y = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  # Add both significance lines
  geom_hline(yintercept = -log10(threshold), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(threshold_2), color = "darkgreen", linetype = "dashed") 


# Results for Nieuwpoort data points
gemma_results_Nieu <- fread("gemma_Np_results.assoc.txt")

gemma_results_Nieu$chr <- as.factor(gemma_results_Nieu$chr) 
gemma_results_Nieu$p_wald <- as.numeric(gemma_results_Nieu$p_wald)

unique(gemma_results_Nieu$chr)

# Create a data frame for chromosome labels and positions
chr_labels <- data.frame(chr = levels(gemma_results_Nieu$chr),
                         label = paste("chr", levels(gemma_results_Nieu$chr), sep = ""))

# Create the Manhattan plot
ggplot(gemma_results_Nieu, aes(x = ps, y = -log10(p_wald))) +
  geom_point(aes(color = chr), size = 1, alpha = 0.6) +
  scale_color_manual(values = rep("blue", length(unique(gemma_results_Nieu$chr)))) +  # Adjust color if needed
  theme_minimal() +
  labs(title = "Manhattan Plot",
       x = "Chromosome Position",
       y = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  # Add both significance lines
  geom_hline(yintercept = -log10(threshold), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(threshold_2), color = "darkgreen", linetype = "dashed")

# Create chromosome mapping based on your provided coordinates
chromosome_mapping <- data.frame(
  chromosome = 1:10,  # Assign chromosome numbers 1 to 10
  scaffold = c("CM008230.1_RagTag", "CM008233.1_RagTag", "CM008234.1_RagTag",
               "CM008235.1_RagTag", "CM008236.1_RagTag", "CM008237.1_RagTag",
               "CM008238.1_RagTag", "CM008239.1_RagTag", "CM008240.1_RagTag",
               "CM008231.1_RagTag"),  
  chromStarts = c(1, 131030212, 222999902, 278345551, 320921871, 
                  379901392, 443804068, 490057698, 517499477, 549277660),
  chromMid = c(61015105.5, 172515056, 246172725.5, 295133710, 345911630.5, 
               407352729, 462430882, 499278586.5, 528888567.5, 561004560.5)
)

gemma_results_Nieu <- gemma_results_Nieu %>%
  mutate(chr = gsub("RagTag", "_RagTag", chr))  # Ensure consistent formatting

chromosome_mapping <- chromosome_mapping %>%
  mutate(scaffold = as.character(scaffold))

# Merge corrected data
gemma_results_Nieu <- gemma_results_Nieu %>%
  left_join(chromosome_mapping, by = c("chr" = "scaffold")) %>%
  mutate(pos_cumulative = ps + chromStarts)

# Check if merge worked
head(gemma_results_Nieu)

library(RColorBrewer)

# Generate a unique color per chromosome
# Color palettes https://www.datanovia.com/en/blog/top-r-color-palettes-to-know-for-great-data-visualization/
chrom_colors <- brewer.pal(n = length(unique(chromosome_mapping$chromosome)), name = "Spectral")

# Create the Manhattan plot
ggplot(gemma_results_Nieu, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(aes(color = factor(chromosome)), size = 1, alpha = 0.6) +
  scale_color_manual(values = chrom_colors) +  # Assign unique colors
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,  # Label at chromosome midpoints
    labels = chromosome_mapping$chromosome  # Show chromosome numbers
  ) +
  theme_minimal() +
  labs(title = "Manhattan Plot, Nieuwpoort data",
       x = "Chromosome",
       y = "-log10(p-value)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Make labels horizontal
        axis.ticks.x = element_blank(),
        legend.position = "none") +
  geom_hline(yintercept = -log10(threshold), color = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(threshold_2), color = "darkgreen", linetype = "dashed")


library(viridis)

# Use better contrast colors
chrom_colors <- viridis(length(unique(chromosome_mapping$chromosome)), option = "D", begin = 0.1, end = 0.9)

ggplot(gemma_results_Nieu, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(aes(color = factor(chromosome), size = -log10(p_wald)), alpha = 0.6) +
  geom_point(data = subset(gemma_results_Nieu, p_wald < threshold),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 0.8) +
  scale_color_manual(values = chrom_colors) +
  scale_size_continuous(range = c(0.5, 2)) +  
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,
    labels = chromosome_mapping$chromosome
  ) +
  theme_minimal() +
  labs(title = "Manhattan Plot: Nieuwpoort Data",
       subtitle = "Genome-wide association results",
       x = "Chromosome",
       y = expression(-log[10](p-value))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic"),
        axis.title = element_text(size = 14)) +
  geom_hline(yintercept = -log10(threshold), color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = -log10(threshold_2), color = "darkgreen", linetype = "dashed", size = 1) +
  annotate("text", x = max(gemma_results_Nieu$pos_cumulative) * 0.95, 
           y = -log10(threshold) + 0.2, label = "Genome-wide threshold", color = "red") +
  annotate("text", x = max(gemma_results_Nieu$pos_cumulative) * 0.95, 
           y = -log10(threshold_2) + 0.2, label = "Suggestive threshold", color = "darkgreen") +
  geom_point(aes(color = factor(chromosome), size = -log10(p_wald)), alpha = 0.6) +
  geom_point(data = subset(gemma_results_Nieu, p_wald < threshold_2),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 0.8)  # Red for significant hits
