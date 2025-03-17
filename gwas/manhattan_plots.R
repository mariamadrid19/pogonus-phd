library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)

# Create chromosome mapping
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

chromosome_mapping <- chromosome_mapping %>%
  mutate(scaffold = as.character(scaffold))

# Results for Heist/Nieuwpoort data points
gemma_results_HeistNpt <- fread("HeistNpt_final.assoc.txt")

gemma_results_HeistNpt$chr <- as.factor(gemma_results_HeistNpt$chr)  # Make sure chromosome is a factor
gemma_results_HeistNpt$p_wald <- as.numeric(gemma_results_HeistNpt$p_wald)
gemma_results_HeistNpt$ps <- as.numeric(gemma_results_HeistNpt$ps)

gemma_results_HeistNpt <- gemma_results_HeistNpt %>%
  left_join(chromosome_mapping, by = c("chr" = "scaffold")) %>%
  mutate(pos_cumulative = ps + chromStarts)

# Apply Bonferroni correction
gemma_results_HeistNpt$p_wald_bonferroni <- p.adjust(gemma_results_HeistNpt$p_wald, method="bonferroni")

# Apply FDR correction (Benjamini-Hochberg)
gemma_results_HeistNpt$p_wald_fdr <- p.adjust(gemma_results_HeistNpt$p_wald, method="fdr")

num_tests <- nrow(gemma_results_HeistNpt)  # Total number of SNPs tested
bonferroni_threshold <- 0.05 / num_tests  # Genome-wide significance threshold
bonferroni_threshold_log <- -log10(bonferroni_threshold)  # Convert to -log10 scale

# New plot to include a Bonferroni significance line
ggplot(gemma_results_HeistNpt, aes(x = pos_cumulative, y = -log10(p_wald))) +
  geom_point(aes(color = factor(chromosome), size = -log10(p_wald)), alpha = 0.6) +
  geom_point(data = subset(gemma_results_HeistNpt, p_wald < bonferroni_threshold),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 0.8) +
  scale_color_manual(values = chrom_colors) +
  scale_size_continuous(range = c(0.5, 2)) +  
  scale_x_continuous(
    breaks = chromosome_mapping$chromMid,
    labels = chromosome_mapping$chromosome
  ) +
  theme_minimal() +
  labs(title = "Manhattan Plot: Heist/Nieuwpoort Data",
       subtitle = "Genome-wide association results",
       x = "Chromosome",
       y = expression(-log[10](p-value))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12, face = "italic"),
        axis.title = element_text(size = 14)) +
  geom_hline(yintercept = bonferroni_threshold_log, color = "blue", linetype = "dashed", size = 0.5) +  # Bonferroni threshold
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed", size = 0.5) +  # Nominal significance
  geom_point(aes(color = factor(chromosome), size = -log10(p_wald)), alpha = 0.6) +
  geom_point(data = subset(gemma_results_HeistNpt, p_wald < bonferroni_threshold),
             aes(x = pos_cumulative, y = -log10(p_wald)),
             color = "red", size = 2, alpha = 0.8) +  # Red for significant hits
  annotate("text", x = max(gemma_results_HeistNpt$pos_cumulative) * 0.95, 
           y = bonferroni_threshold_log + 0.2, label = "Bonferroni Threshold", color = "blue") +
  annotate("text", x = max(gemma_results_HeistNpt$pos_cumulative) * 0.95, 
           y = -log10(0.05) + 0.2, label = "p < 0.05 Threshold", color = "red")
