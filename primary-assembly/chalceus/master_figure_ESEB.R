library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales) 
library(colorspace)
library(patchwork)
library(ggfortify)
library(gridExtra)
library(ggpubr) 
library(scales)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

setwd("~/Documents/phd-pogonus/master-figure")

#define maximum LG
LGmax=20

chrom_data <- read.table("scaflengths.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
features <- read.table("map8_ann.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
map <- read.table("map8_all_ann.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

chrom_data <- chrom_data %>%  mutate(start = lag(cumsum(length), default = 0),end = start + length,color = rep(c("gray80", "gray40"), length.out = n()))
chrom_data_chr <- chrom_data %>% filter(grepl("^CHR", chrom))

# Convert category to a factor
features <- subset(features,LG<LGmax)
features$LG <- as.factor(features$LG)

# Merge to get chromosome start positions from `chrom_data_chr`
features_merged <- merge(
  features, 
  chrom_data_chr[, c("chrom", "start")],
  by.x = "chr", 
  by.y = "chrom"
)
map_merged <- merge(map, chrom_data_chr[, c("chrom", "start")], by = "chrom")

# Compute global genome coordinate
features_merged$genome_pos <- as.numeric(features_merged$start) + as.numeric(features_merged$snp)
map_merged$genome_pos <- as.numeric(map_merged$start) + as.numeric(map_merged$pos)
map_merged$cM_f <- as.numeric(map_merged$cM_f)
map_merged$cM_m <- as.numeric(map_merged$cM_m)
str(map_merged)

# Plot LG
plot = ggplot() + geom_rect(data = chrom_data_chr, aes(xmin = start, xmax = end, ymin = 0, ymax = LGmax, fill = color),alpha=0.2)+scale_fill_identity()
plot = plot+geom_point(data = features_merged, mapping=aes(x = genome_pos, y = LG), shape=21, size = 1, fill="black", alpha = 0.2)
plot = plot+theme_minimal()
plot

# Plot map
max_cM<-max(map_merged$cM_f)
plot = ggplot() + geom_rect(data = chrom_data_chr, aes(xmin = start, xmax = end, ymin = 0, ymax = max_cM, fill = color),alpha=0.2)+scale_fill_identity()
plot = plot+geom_point(data = map_merged, mapping=aes(x = genome_pos, y = cM_f,color=LG), size = 1, alpha = 0.5)
plot = plot+theme_minimal()+xlim(0,500000000)
plot

# Define max cM again in case it's changed
max_cM <- max(map_merged$cM_f)
map_merged_no_LG12 <- subset(map_merged, LG != "LG12")

# Plot
plot <- ggplot() +
  geom_rect(data = chrom_data_chr,
            aes(xmin = start, xmax = end, ymin = 0, ymax = max_cM, fill = color),
            alpha = 0.2) +
  scale_fill_identity() +
  geom_point(data = map_merged_no_LG12,
             aes(x = genome_pos, y = cM_f, color = LG),
             size = 1, alpha = 0.5) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = (chrom_data_chr$start + chrom_data_chr$end) / 2,
    labels = 1:nrow(chrom_data_chr),
    limits = c(min(chrom_data_chr$start), max(chrom_data_chr$end))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(y = "Genetic Position (cM)")
plot


#define maximum LG
LGmax=20

chrom_data <- read.table("scaflengths.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
features <- read.table("map8_ann.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
map <- read.table("map8_all_ann_mrecom0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

chrom_data <- chrom_data %>%  mutate(start = lag(cumsum(length), default = 0),end = start + length,color = rep(c("gray80", "gray40"), length.out = n()))
chrom_data_chr <- chrom_data %>% filter(grepl("^CHR", chrom))

# Convert category to a factor
features <- subset(features,LG<LGmax)
features$LG <- as.factor(features$LG)

# Merge to get chromosome start positions from `chrom_data_chr`
features_merged <- merge(
  features, 
  chrom_data_chr[, c("chrom", "start")],
  by.x = "chr", 
  by.y = "chrom"
)
map_merged <- merge(map, chrom_data_chr[, c("chrom", "start")], by = "chrom")

# Compute global genome coordinate
features_merged$genome_pos <- as.numeric(features_merged$start) + as.numeric(features_merged$snp)
map_merged$genome_pos <- as.numeric(map_merged$start) + as.numeric(map_merged$pos)
map_merged$cM_f <- as.numeric(map_merged$cM_f)
map_merged$cM_m <- as.numeric(map_merged$cM_m)
str(map_merged)

# Plot
plot <- ggplot() +
  geom_rect(data = chrom_data_chr,
            aes(xmin = start, xmax = end, ymin = 0, ymax = max_cM, fill = color),
            alpha = 0.2) +
  scale_fill_identity() +
  geom_point(data = map_merged_no_LG12,
             aes(x = genome_pos, y = cM_f, color = LG),
             size = 1, alpha = 0.5) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = (chrom_data_chr$start + chrom_data_chr$end) / 2,
    labels = 1:nrow(chrom_data_chr),
    limits = c(min(chrom_data_chr$start), max(chrom_data_chr$end))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  labs(y = "Genetic Position (cM)")
plot


spain_stats <- read.table("spain_stats.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
scaf_coords <- read.table("scaf_coords.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
chrom_colors <- rep(c("lightgray", "darkgray"), length.out = 11)

# Merge and drop duplicate columns
genome_data <- merge(
  chrom_data_chr,
  scaf_coords,
  by.x = "chrom",
  by.y = "scaffold"
)

# Remove unwanted duplicate columns
genome_data <- subset(genome_data, select = -c(length.y, chrom_id, chromLengths))

# Rename length.x
names(genome_data)[names(genome_data) == "length.x"] <- "length"

# Sort by chromosome number
genome_data$chrom_num <- as.numeric(sub("CHR", "", genome_data$chrom))
genome_data <- genome_data[order(genome_data$chrom_num), ]
genome_data$chrom_num <- NULL

genome_data
#write.table(genome_data, file = "genome_data.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


# Fst plot, Spain FST
Fst <- ggplot(spain_stats, aes(x = chromPos, y = FstWC_S_HUE_T_S_COT_S_egg)) +
  geom_point(aes(color = factor(chromosome))) +
  scale_color_manual(values = chrom_colors) +
  scale_x_continuous(
    breaks = genome_data$chromMid,
    labels = genome_data$chromosome
  ) +
  scale_y_continuous(
    breaks = function(x) {
      brks <- pretty(x)
      brks[seq(1, length(brks), by = 2)]
    }
  ) +
  theme_minimal() +
  labs(x = "Chromosome",
       y = expression(F[st])) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16))

Fst


# Plot linkage map with new genome_data dataset
linkage_map <- ggplot() +
  geom_rect(data = genome_data,
            aes(xmin = start, xmax = end, ymin = 0, ymax = max_cM, fill = color),
            alpha = 0.2) +
  scale_fill_identity() +
  geom_point(data = map_merged_no_LG12,
             aes(x = genome_pos, y = cM_f, color = LG),
             size = 1, alpha = 0.5) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = (genome_data$start + genome_data$end) / 2,
    labels = 1:nrow(genome_data),
    limits = c(min(genome_data$start), max(genome_data$end))
  ) +
  theme_minimal()  +
  labs(x = "Chromosome",
       y = "Genetic position (cM)") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16))

linkage_map

# Remove x-axis title and labels from top plot
linkage_map_clean <- linkage_map +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

# Stack plots with patchwork
combined_plot <- linkage_map_clean / Fst  +
  plot_layout(heights = c(1, 1)) # adjust proportions if needed

combined_plot



# COVERAGE DATA
# Create list to store all processed data frames
setwd("~/Documents/phd-pogonus/master-figure/raw-cov")
# Define sample-to-population mapping
sample_info <- data.frame(
  sample = c(
    "GC129388", "GC129389", "GC129390", "GC129391", "GC129392", "GC129393",
    "GC129394", "GC129395", "GC129396", "GC129397", "GC129398", "GC129399",
    "GC129400", "GC129401", "GC129402", "GC129403", "GC129404", "GC129405",
    "GC129406", "GC129407", "GC129408", "GC129409", "GC129410", "GC129411",
    "GC129412", "GC129413", "GC129414", "GC129415", "GC129416", "GC129417",
    "GC129418", "GC129419", "GC129420", "GC129421", "GC129422", "GC129423",
    "GC129424", "GC129425", "GC129426", "GC136078", "GC136079", "GC136080",
    "GC136081", "GC136082", "GC136083", "GC136084", "GC136085", "GC136086",
    "GC136087", "GC136088", "GC136089", "GC136090", "GC136091", "GC136092",
    "GC136093", "GC136094", "GC136095", "GC136096", "GC136097", "GC136098",
    "GC136099", "GC136100", "GC136101", "GC136102", "GC136103", "GC136104",
    "GC136105", "GC136106", "GC136107", "GC136108", "GC136109", "GC136110",
    "GC136111", "GC136112", "GC136113", "GC136114", "GC136115", "GC136116"
  ),
  population = c(
    "B_NIE_T", "B_NIE_T", "B_NIE_T", "B_NIE_T", "B_NIE_T", "B_NIE_T",
    "B_DUD_S", "B_DUD_S", "B_DUD_S", "B_DUD_S", "B_DUD_S", "B_DUD_S",
    "F_GUE_T", "F_GUE_T", "F_GUE_T", "F_GUE_T", "F_GUE_T", "F_GUE_T",
    "F_GUE_S", "F_GUE_S", "F_GUE_S", "F_GUE_S", "F_GUE_S", "F_GUE_S",
    "P_AVE_T", "P_AVE_T", "P_AVE_T", "P_AVE_T", "P_AVE_T", "P_AVE_S",
    "P_AVE_S", "P_AVE_S", "P_AVE_S", "P_AVE_S", "F_CAM_S", "F_CAM_S",
    "E_SEV_T", "E_SEV_T", "E_SEV_T", "F_GUE_T", "F_GUE_T", "F_GUE_T",
    "F_GUE_T", "F_GUE_T", "F_GUE_T", "B_DUD_S", "B_DUD_S", "B_DUD_S",
    "B_DUD_S", "B_DUD_S", "B_DUD_S", "B_NIE_T", "B_NIE_T", "B_NIE_T",
    "B_NIE_T", "B_NIE_T", "B_NIE_T", "F_GUE_S", "F_GUE_S", "F_GUE_S",
    "F_GUE_S", "F_GUE_S", "F_GUE_S", "F_CAM_S", "F_CAM_S", "F_CAM_S",
    "E_SEV_T", "E_SEV_T", "S_HUE_T", "S_HUE_T", "S_COT_S", "S_COT_S",
    "S_COT_S", "S_COT_S", "S_COT_S", "S_HUE_T", "S_HUE_T", "S_HUE_T"
  ),
  stringsAsFactors = FALSE
)

sample_info$ecotype <- ifelse(
  grepl("_T$", sample_info$population), "T",
  ifelse(grepl("_S$", sample_info$population), "S", NA)
)

# Create list to store all processed data frames
processed_list <- list()

# Loop over each sample
for (i in 1:nrow(sample_info)) {
  sample_id <- sample_info$sample[i]
  pop <- sample_info$population[i]
  eco <- sample_info$ecotype[i]
  
  # Read the .raw.pcov file
  file_path <- paste0(sample_id, ".raw.pcov")
  data <- read.delim(file_path, header = FALSE, comment.char = "#")
  colnames(data) <- c("region", "coverage", "length")
  
  # Process the data
  data <- data %>%
    filter(grepl(">", region)) %>%
    mutate(region = gsub(">", "", region)) %>%
    filter(str_starts(region, "CHR")) %>%
    mutate(
      chromosome = as.numeric(str_extract(region, "\\d+")),  # numeric chr
      region = gsub(".*_", "", region),
      region = as.numeric(region),
      length = as.numeric(length),
      position = region * length
    ) %>%
    arrange(chromosome) %>%
    group_by(chromosome) %>%
    ungroup()
  
  # Merge with genome_data instead of chrom_coords
  merged_df <- merge(
    data,
    genome_data[, c("chromosome", "chromStarts")],  # only needed cols
    by = "chromosome",
    all.x = TRUE
  )
  
  # Rename columns
  names(merged_df)[names(merged_df) == "length"] <- "window_length"
  
  # Compute genomic position and add metadata
  merged_df$genomic_position <- merged_df$position + merged_df$chromStarts - 1
  merged_df$individual <- sample_id
  merged_df$population <- pop
  merged_df$ecotype <- eco
  average_coverage <- mean(merged_df$coverage, na.rm = TRUE)
  merged_df$normalized_coverage <- merged_df$coverage / average_coverage
  
  # Save in global env and list
  assign(sample_id, merged_df, envir = .GlobalEnv)
  processed_list[[sample_id]] <- merged_df
}

# Merge all processed sample data frames
merged_data <- do.call(rbind, processed_list)

# Create an empty list to hold plots
plots <- list()

# Get all chromosomes present in the merged data
chromosomes <- sort(unique(merged_data$chromosome))

# Loop through chromosomes and make a plot for each
for (chr in chromosomes) {
  chr_data <- merged_data %>% filter(chromosome == chr)
  
  p <- ggplot(chr_data, aes(
    x = genomic_position / 1e6,  # convert to Mbp
    y = normalized_coverage,
    group = individual
  )) +
    geom_smooth(
      aes(color = ecotype),
      method = "loess",
      span = 0.2,   # smoothness
      se = FALSE,
      size = 0.7
    ) +
    scale_color_manual(
      values = c("T" = alpha("blue", 0.3), "S" = alpha("red", 0.3))
    ) +
    coord_cartesian(ylim = c(0, 3)) +
    theme_minimal() +
    labs(
      title = paste("Chromosome", chr),
      x = "Genomic position (Mbp)",
      y = "Coverage"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold")
    )
  
  plots[[length(plots) + 1]] <- p
}
# Arrange all chromosome plots in a grid
grid.arrange(grobs = plots, ncol = 4)



# Ensure ecotype is a factor for consistent coloring
merged_data <- merged_data %>%
  mutate(ecotype = factor(ecotype, levels = c("T", "S")))

# Create an empty list to hold plots
plots <- list()

# Get all chromosomes present in the merged data
chromosomes <- sort(unique(merged_data$chromosome))

# Coverage plot with chromosome-aligned backgrounds
coverage_plot <- ggplot() +
  # Grey chromosome backgrounds
  geom_rect(
    data = genome_data,
    aes(
      xmin = chromStarts, xmax = chromEnds,
      ymin = 0, ymax = 3,  # matches coord_cartesian y limits
      fill = color
    ),
    alpha = 0.2
  ) +
  scale_fill_identity() +
  
  # Smoothed coverage lines per individual per chromosome
  geom_smooth(
    data = merged_data,
    aes(
      x = genomic_position,
      y = normalized_coverage,
      group = interaction(individual, chromosome),
      color = ecotype
    ),
    method = "loess",
    span = 0.2,
    se = FALSE,
    size = 0.5
  ) +
  
  # Colors with transparency for ecotypes
  scale_color_manual(
    values = c("T" = alpha("blue", 0.3), "S" = alpha("red", 0.3))
  ) +
  
  # Same X axis as Fst & linkage_map
  scale_x_continuous(
    name = "Chromosome",
    breaks = genome_data$chromMid,
    labels = genome_data$chromosome,
    limits = c(min(genome_data$chromStarts), max(genome_data$chromEnds))
  ) +
  
  coord_cartesian(ylim = c(0, 2)) +
  theme_minimal() +
  labs(
    y = "Normalized coverage"
  ) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = c(0.09, 0.91),
    legend.direction = "horizontal",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16)
  )

coverage_plot_clean <- coverage_plot +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
coverage_plot_clean

# Stack plots with patchwork
combined_plot <- linkage_map_clean  / coverage_plot_clean / Fst +
  plot_layout(heights = c(1, 1, 1)) # adjust proportions if needed

combined_plot


# plot scaffolds
setwd("~/Documents/phd-pogonus/master-figure")
agp <- read.table('Pchalceus_SW.sorted.fasta.agp')
agp <- subset(agp, agp$V5 == 'W')

rows_with_pattern <- apply(agp, 1, function(row) any(grepl('CHR', row)))
agp <- agp[rows_with_pattern, ]

colnames(agp) <- c('chrom', 'start', 'end', 'part_number','component_type','component_id','component_beg','component_end','orientation')

# Ensure chrom is character in both datasets
agp$chrom <- as.character(agp$chrom)
genome_data$chrom <- as.character(genome_data$chrom)

# Merge AGP with genome_data
agpM <- merge(
  agp,
  genome_data[, c("chrom", "chromStarts")],
  by = "chrom",
  all.x = TRUE
)

# Compute genome-wide start and end positions
agpM$chromPosStart <- agpM$start + agpM$chromStarts - 1
agpM$chromPosEnd   <- agpM$end   + agpM$chromStarts - 1
agpM

# Create alternating scaffold index within each chromosome
agpM <- agpM %>%
  arrange(chrom, chromPosStart) %>%
  group_by(chrom) %>%
  mutate(scaffold_index = row_number() %% 2) %>%
  ungroup()

# --- FST PLOT with gradient and scale bar ---
Fst_plot <- ggplot(spain_stats, aes(x = chromPos, y = FstWC_S_HUE_T_S_COT_S_egg)) +
  geom_point(aes(color = FstWC_S_HUE_T_S_COT_S_egg)) +  # color by Fst value
  scale_color_gradient(low = "black", high = "red", name = expression(F[st])) + # gradient legend
  scale_x_continuous(
    breaks = genome_data$chromMid,
    labels = genome_data$chromosome,
    limits = c(min(genome_data$chromStarts), max(genome_data$chromEnds))
  ) +
  scale_y_continuous(
    breaks = function(x) {
      brks <- pretty(x)
      brks[seq(1, length(brks), by = 2)]
    }
  ) +
  theme_minimal() +
  labs(x = "Chromosome",
       y = expression(F[st])) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16)) +
  # Add scale bar (20 Mb)
  annotate("segment",
           x = 30e6, xend = 50e6,  # position in genome coordinates (50 Mb → 70 Mb)
           y = 0.9, yend = 0.9,    # vertical position
           colour = "black", size = 1.2) +
  annotate("text",
           x = 40e6, y = 0.85, label = "20 Mb",
           size = 4, hjust = 0.5)

Fst_plot

# --- SCAFFOLD TRACK PLOT ---
scaffold_plot <- ggplot() +
  geom_rect(
    data = agpM,
    aes(
      xmin = chromPosStart / 1e6,
      xmax = chromPosEnd / 1e6,
      ymin = 0,
      ymax = 1,
      fill = factor(scaffold_index)
    ),
    color = NA
  ) +
  scale_fill_manual(values = c("darkturquoise", "deepskyblue4")) +
  scale_x_continuous(
    breaks = genome_data$chromMid / 1e6,
    labels = genome_data$chromosome
  ) +
  theme_minimal() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )
scaffold_plot
    
# --- COMBINE EVERYTHING
combined_plot <- linkage_map_clean  / scaffold_plot / coverage_plot_clean / Fst_plot +
  plot_layout(heights = c(1, 0.1, 1, 1))

combined_plot

combined_plot_2 <- linkage_map_clean  / scaffold_plot / coverage_plot_clean / Fst +
  plot_layout(heights = c(1, 0.1, 1, 1))
combined_plot_2


gemma_results_NP25_1 <- read.table("gemma_results_NP25_MRWS.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
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
  theme_minimal() +
  labs(title = "%MRWS (n=192)",
       x = NULL,
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_pwald_1

gemma_results_NP25_3 <- read.table("gemma_results_NP25_chr5AA.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

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
  theme_minimal() +
  labs(title = "%MRWS - chr5:AA individuals (n=120)",
       x = "Chromosome",
       y = expression(-log[10](p))) +
  theme(axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_text(size = 16))
manhattan_pwald_3

Fst_plot_clean <- ggplot(spain_stats, aes(x = chromPos, y = FstWC_S_HUE_T_S_COT_S_egg)) +
  geom_point(aes(color = FstWC_S_HUE_T_S_COT_S_egg)) +  # color by Fst value
  scale_color_gradient(low = "black", high = "red", name = expression(F[st])) + # gradient legend
  scale_x_continuous(
    breaks = genome_data$chromMid,
    labels = genome_data$chromosome,
    limits = c(min(genome_data$chromStarts), max(genome_data$chromEnds))
  ) +
  scale_y_continuous(
    breaks = function(x) {
      brks <- pretty(x)
      brks[seq(1, length(brks), by = 2)]
    }
  ) +
  theme_minimal() +
  labs(x = NULL,
       y = expression(F[st])) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 16)) +
  # Add scale bar (20 Mb)
  annotate("segment",
           x = 30e6, xend = 50e6,  # position in genome coordinates (50 Mb → 70 Mb)
           y = 0.9, yend = 0.9,    # vertical position
           colour = "black", size = 1.2) +
  annotate("text",
           x = 40e6, y = 0.85, label = "20 Mb",
           size = 4, hjust = 0.5)

Fst_plot_clean

combined_plot_GWAS <- Fst_plot_clean/ manhattan_pwald_1 / manhattan_pwald_3 +
  plot_layout(heights = c(0.5, 1, 1))
combined_plot_GWAS


### PLOT GENE DENSITY BASED ON EVIANN ANNOTATION FILE 

exon_density <- read.table("exon_density_100kb.tsv", header = FALSE, sep = "\t")
colnames(exon_density) <- c("chrom", "start", "end", "exon_percent")
exon_density$mid <- (exon_density$start + exon_density$end) / 2

# Align exon_density positions to genome-wide coordinates
exon_density_aligned <- exon_density %>%
  left_join(genome_data %>% select(chrom, chromStarts), by = "chrom") %>%
  mutate(
    genomic_position = mid + chromStarts
  )

# Exon density plot with chromosome breaks + smoothing
exon_density_plot <- ggplot() +
  # Grey chromosome backgrounds
  geom_rect(
    data = genome_data,
    aes(
      xmin = chromStarts, xmax = chromEnds,
      ymin = 0, ymax = 30,   # fixed y limit
      fill = color
    ),
    alpha = 0.2
  ) +
  scale_fill_identity() +
  
  # Smoothed exon density per chromosome
  geom_smooth(
    data = exon_density_aligned,
    aes(
      x = genomic_position,
      y = exon_percent,
      group = chrom
    ),
    method = "loess",
    span = 0.1,
    se = FALSE,
    size = 1,
    color = "darkgreen"
  ) +
  
  # Same X axis as other genome-wide plots
  scale_x_continuous(
    name = NULL,
    breaks = genome_data$chromMid,
    labels = genome_data$chromosome,
    limits = c(min(genome_data$chromStarts), max(genome_data$chromEnds))
  ) +
  
  coord_cartesian(ylim = c(0, 30)) +
  
  theme_minimal() +
  labs(
    y = "Gene %"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16)
  )
exon_density_plot

combined_plot_genes <- Fst_plot_clean / exon_density_plot +
  plot_layout(heights = c(1, 1))
combined_plot_genes

# Step 1 — Mapping table
lg_to_chr <- data.frame(
  LG = c("LG01","LG02","LG03","LG04","LG05","LG06","LG07","LG08","LG09","LG10","LG11"),
  chrom = c("CHR1","CHR2","CHR4","CHR5","CHR3","CHR6","CHR7","CHR8","CHR9","CHR10","CHR11")
)

# Step 2: Ensure no leftover chrom column
map_merged_no_LG12_fixed <- map_merged_no_LG12 %>%
  select(-any_of("chrom")) %>%  # Remove if it exists
  left_join(lg_to_chr, by = "LG") %>%  # Add correct chrom
  left_join(genome_data %>% select(chrom, chromStarts), by = "chrom") %>%
  mutate(
    pos = as.numeric(pos),
    chromStarts = as.numeric(chromStarts),
    genome_pos_aligned = pos + chromStarts,
    chrom = factor(chrom, levels = paste0("CHR", 1:11))  # Force CHR order
  )

# Step 3: Plot using the fixed data
linkage_map_clean <- ggplot() +
  geom_rect(
    data = genome_data,
    aes(
      xmin = chromStarts,
      xmax = chromEnds,
      ymin = 0,
      ymax = max_cM,
      fill = color
    ),
    alpha = 0.2
  ) +
  scale_fill_identity() +
  geom_point(
    data = map_merged_no_LG12_fixed,
    aes(
      x = genome_pos_aligned,
      y = cM_f
    ),
    color = "gray30",
    size = 0.6,
    alpha = 0.5
  ) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = genome_data$chromMid,
    labels = genome_data$chromosome,
    limits = c(min(genome_data$chromStarts), max(genome_data$chromEnds))
  ) +
  theme_minimal() +
  labs(y = "Genetic position (cM)") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16)
  )

linkage_map_clean

Fst_plot <- ggplot() +
  geom_rect(
    data = genome_data,
    aes(
      xmin = chromStarts,
      xmax = chromEnds,
      ymin = 0,
      ymax = max(spain_stats$FstWC_S_HUE_T_S_COT_S_egg, na.rm = TRUE),
      fill = color
    ),
    alpha = 0.2
  ) +
  scale_fill_identity() +
  # Fst points
  geom_point(
    data = spain_stats,
    aes(x = chromPos, y = FstWC_S_HUE_T_S_COT_S_egg, color = FstWC_S_HUE_T_S_COT_S_egg)
  ) +
  scale_color_gradient(low = "black", high = "red", name = expression(F[st])) +
  # X axis
  scale_x_continuous(
    breaks = genome_data$chromMid,
    labels = genome_data$chromosome,
    limits = c(min(genome_data$chromStarts), max(genome_data$chromEnds))
  ) +
  # Y axis
  scale_y_continuous(
    breaks = function(x) {
      brks <- pretty(x)
      brks[seq(1, length(brks), by = 2)]
    }
  ) +
  theme_minimal() +
  labs(x = "Chromosome", y = expression(F[st])) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16)
  ) +
  # Add scale bar (20 Mb)
  annotate(
    "segment",
    x = 30e6, xend = 50e6,
    y = 0.9, yend = 0.9,
    colour = "black", size = 1.2
  ) +
  annotate(
    "text",
    x = 40e6, y = 0.85, label = "20 Mb",
    size = 4, hjust = 0.5
  )
Fst_plot

combined_plot_3 <- linkage_map_clean  / scaffold_plot / coverage_plot_clean / exon_density_plot / Fst_plot +
  plot_layout(heights = c(0.8, 0.1, 0.6, 1, 1))
combined_plot_3

setwd("~/Documents/phd-pogonus/master-figure")

# Read the GFF file into a data frame
gff_file <- "Pchalceus_SW.sorted.fasta.out.gff"

gff_data <- read.table(gff_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)
colnames(gff_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gff_data <- gff_data %>%
  filter(grepl("CHR", seqname))


# Create GRanges object
gr <- GRanges(seqnames = gff_data$seqname,
              ranges = IRanges(start = gff_data$start, end = gff_data$end),
              strand = gff_data$strand)

# Get unique scaffolds
scaffolds <- unique(seqnames(gr))
scaffolds <- unique(scaffolds[grep("CHR", scaffolds)])[c(1:11)]

result_list <- lapply(scaffolds, function(scaffold) {
  # Subset the GRanges object for the current scaffold
  gr_scaffold <- gr[seqnames(gr) == scaffold]
  
  # Determine the length of the scaffold
  scaffold_length <- max(end(gr_scaffold))
  
  # Create sliding windows for the scaffold
  windows <- slidingWindows(GRanges(seqnames = scaffold, ranges = IRanges(start = 1, end = scaffold_length)), 
                            width = 100000, step = 100000)[[1]]
  
  # Calculate the amount covered by features in each window
  coverage <- sapply(seq_along(windows), function(i) {
    window <- windows[i]
    overlaps <- pintersect(rep(window, length(gr_scaffold)), gr_scaffold)
    sum(width(overlaps))
  })
  
  # Create a data frame for the results
  data.frame(
    seqnames = scaffold,
    start = start(windows),
    end = end(windows),
    coverage = coverage
  )
})

# Combine the results for all scaffolds
result <- do.call(rbind, result_list)

# Print the result
print(result)

result$mid <- (result$start+result$end)/2
result$perc <- result$coverage/100000

scaffolds <- c("CHR1", "CHR2", "CHR3", "CHR4", "CHR5", "CHR6", "CHR7", "CHR8","CHR9","CHR10","CHR11")

for (scaffold in scaffolds) {
  result_scaffolds <- subset(result, result$seqnames == scaffold)
  
  plot(result_scaffolds$start, result_scaffolds$perc, pch=19, cex = 0.3, col="#CD9600",main = scaffold, xlab = "Chromosome Position",ylab="Repeat Content")
}


#add to combined plot
merged_result <- result %>%
  left_join(genome_data, by = c("seqnames" = "chrom"))

merged_result$perc <- pmin(merged_result$perc, 1)

merged_result <- merged_result %>%
  mutate(genomic_position = mid + chromStarts)

merged_result$perc <- pmin(merged_result$perc, 1) * 100

repeat_plot <- ggplot() +
  geom_rect(
    data = genome_data,
    aes(
      xmin = chromStarts,
      xmax = chromEnds,
      ymin = 0,
      ymax = 100,
      fill = color
    ),
    alpha = 0.2,
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_identity() +
  geom_point(
    data = merged_result,
    aes(x = genomic_position, y = perc),
    color = "#CD9600",
    size = 0.5
  ) +
  scale_x_continuous(
    name = "Chromosome",
    breaks = genome_data$chromMid,
    labels = genome_data$chromosome,
    limits = c(min(genome_data$chromStarts), max(genome_data$chromEnds))
  ) +
  scale_y_continuous(
    name = "Repeat %",
    limits = c(0, 100),
    breaks = c(0, 50, 100)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 16)
  )
repeat_plot_clean <- repeat_plot +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )
repeat_plot_clean


combined_plot_4 <- linkage_map_clean  / scaffold_plot / coverage_plot_clean / repeat_plot_clean / exon_density_plot / Fst_plot +
  plot_layout(heights = c(0.8, 0.1, 0.6, 0.4, 0.4, 1))
combined_plot_4

#write.table(merged_result, file = "repeat_content.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
