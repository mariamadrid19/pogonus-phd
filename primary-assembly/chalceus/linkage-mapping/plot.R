library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(scales) 
library(colorspace)

#define maximum LG
LGmax=20

chrom_data <- read.table("scaflengths_broken.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
features <- read.table("map8_js4_ann.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
map <- read.table("map8_all_ann.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

chrom_data <- chrom_data %>%  mutate(start = lag(cumsum(length), default = 0),end = start + length,color = rep(c("gray80", "gray40"), length.out = n()))

# Convert category to a factor
features <- subset(features,LG<LGmax)
features$LG <- as.factor(features$LG)

# Merge to get chromosome start positions from `chrom_data`
features_merged <- merge(features, chrom_data[, c("chrom", "start")], by = "chrom")
map_merged <- merge(map, chrom_data[, c("chrom", "start")], by = "chrom")

# Compute global genome coordinate
features_merged$genome_pos <- as.numeric(features_merged$start) + as.numeric(features_merged$pos)
map_merged$genome_pos <- as.numeric(map_merged$start) + as.numeric(map_merged$pos)

# Only retain LGs that are below LGmax


# Plot LG
plot = ggplot() + geom_rect(data = chrom_data, aes(xmin = start, xmax = end, ymin = 0, ymax = LGmax, fill = color),alpha=0.2)+scale_fill_identity()
plot = plot+geom_point(data = features_merged, mapping=aes(x = genome_pos, y = LG), shape=21, size = 1, fill="black", alpha = 0.2)
plot = plot+theme_minimal()
plot

# Plot map
max_cM<-max(map_merged$cM_f)
plot = ggplot() + geom_rect(data = chrom_data, aes(xmin = start, xmax = end, ymin = 0, ymax = max_cM, fill = color),alpha=0.2)+scale_fill_identity()
plot = plot+geom_point(data = map_merged, mapping=aes(x = genome_pos, y = cM_f,color=LG), size = 1, alpha = 0.5)
plot = plot+theme_minimal()+xlim(0,500000000)
plot

# Limit to scaffold_1 to scaffold_20
scaffolds_to_keep <- paste0("scaffold_", 1:20)
chrom_data_filtered <- chrom_data %>%
  filter(chrom %in% scaffolds_to_keep) %>%
  mutate(mid = start + (length / 2))

map_merged_filtered <- map_merged %>%
  filter(chrom %in% scaffolds_to_keep)

# Define max cM again in case it's changed
max_cM <- max(map_merged_filtered$cM_f)

# Plot
plot <- ggplot() +
  geom_rect(data = chrom_data_filtered,
            aes(xmin = start, xmax = end, ymin = 0, ymax = max_cM, fill = color),
            alpha = 0.2) +
  scale_fill_identity() +
  geom_point(data = map_merged_filtered,
             aes(x = genome_pos, y = cM_f, color = LG),
             size = 1, alpha = 0.5) +
  scale_x_continuous(
    name = "Scaffold",
    breaks = chrom_data_filtered$mid,
    labels = chrom_data_filtered$chrom,
    limits = c(min(chrom_data_filtered$start), max(chrom_data_filtered$end))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  labs(y = "Genetic Position (cM)")

plot


# === 1. Read and parse the file ===
map <- read.table("map.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(map) <- c("scaffold", "scaffold_range", "lg_range", "strand", "LG")

# === 2. Parse coordinate ranges ===
parse_coords <- function(coord) {
  as.numeric(str_split(coord, "-", simplify = TRUE))
}

# Scaffold start and end
scaffold_coords <- t(sapply(map$scaffold_range, parse_coords))
map$scaffold_start <- scaffold_coords[,1]
map$scaffold_end <- scaffold_coords[,2]

# Use scaffold midpoint for plotting
map$scaffold_mid <- (map$scaffold_start + map$scaffold_end) / 2

# Orientation
map$orientation <- ifelse(grepl("\\*$", map$lg_range), "reverse", "forward")

# === 3. Prepare plot data ===
# Normalize scaffold positions within LGs
map <- map %>%
  group_by(LG) %>%
  arrange(as.numeric(str_extract(lg_range, "^\\d+")), .by_group = TRUE) %>%
  mutate(xmin = cumsum(lag(scaffold_end - scaffold_start + 1, default = 0)),
         xmax = xmin + (scaffold_end - scaffold_start + 1),
         scaffold_label = str_remove(scaffold, "^scaffold_")) %>%
  ungroup()

# === 4. Plot ===
ggplot(map, aes(xmin = xmin / 1e6, xmax = xmax / 1e6, ymin = LG - 0.4, ymax = LG + 0.4, fill = orientation)) +
  geom_rect(color = "black") +
  scale_y_continuous(breaks = unique(map$LG)) +
  scale_fill_manual(values = c("forward" = "#1f77b4", "reverse" = "#ff7f0e")) +
  labs(
    x = "Position within LG (Mb)",
    y = "Linkage Group",
    title = "Scaffold arrangement by linkage group",
    fill = "Orientation"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )


# === 2. Parse coordinate ranges ===
parse_coords <- function(coord) {
  as.numeric(str_split(coord, "-", simplify = TRUE))
}

# Extract scaffold start/end
scaffold_coords <- t(sapply(map$scaffold_range, parse_coords))
map$scaffold_start <- scaffold_coords[,1]
map$scaffold_end <- scaffold_coords[,2]

# Compute width and midpoints (optional)
map$scaffold_width <- map$scaffold_end - map$scaffold_start + 1

# === 3. Assign scaffold segments to relative x positions within each LG ===
map <- map %>%
  group_by(LG) %>%
  arrange(scaffold, .by_group = TRUE) %>%
  mutate(
    xmin = cumsum(lag(scaffold_width, default = 0)),
    xmax = xmin + scaffold_width,
    scaffold_label = str_remove(scaffold, "^scaffold_")
  ) %>%
  ungroup()

# === 4. Plot, coloring by scaffold ===
library(Polychrome)

# Generate 23 distinct colors with high visual contrast
scaffold_names <- unique(map$scaffold)
scaffold_colors <- createPalette(length(scaffold_names), seedcolors = c("#FF0000", "#00FF00", "#0000FF"))

names(scaffold_colors) <- scaffold_names

# Plot
ggplot(map, aes(xmin = xmin / 1e6, xmax = xmax / 1e6, ymin = LG - 0.4, ymax = LG + 0.4, fill = scaffold)) +
  geom_rect(color = "black") +
  scale_y_continuous(breaks = sort(unique(map$LG))) +
  scale_fill_manual(values = scaffold_colors, name = "Scaffold") +
  labs(
    x = "Position within Linkage Group (Mb)",
    y = "Linkage Group",
    title = "Scaffold Composition of Linkage Groups"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 9)
  )




## FINAL AGP
# Load AGP file
agp <- read.table("final.fasta.agp", header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Name columns according to AGP format
colnames(agp) <- c("chr", "start", "end", "part", "type", "scaffold", "scaf_start", "scaf_end", "orientation")

# Separate components (W) and gaps (U)
components <- agp %>% filter(type == "W")
gaps <- agp %>% filter(type == "U")

# Extract scaffold number for labeling
components$scaffold_num <- str_replace(components$scaffold, "scaffold_", "")

# Assign colors for chromosomes
chrom_colors <- scales::hue_pal()(length(unique(components$chr)))
names(chrom_colors) <- unique(components$chr)

# Add fill color based on chromosome and orientation
components <- components %>%
  mutate(
    chr_color = chrom_colors[chr],
    fill = ifelse(orientation == "+", chr_color, scales::muted(chr_color))  # darker shade for reverse
  )

# Gaps are light gray
gaps$fill <- "lightgray"

# Combine for plotting
plot_data <- bind_rows(
  components %>% select(chr, start, end, fill, scaffold_num),
  gaps %>% mutate(scaffold_num = NA) %>% select(chr, start, end, fill, scaffold_num)
)

# Order chromosomes
plot_data$chr <- factor(plot_data$chr, levels = paste0("CHR", 1:12))

ggplot(plot_data) +
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6,
                ymin = as.numeric(chr) - 0.4, ymax = as.numeric(chr) + 0.4,
                fill = fill), color = "black") +
  geom_text(
    data = filter(plot_data, !is.na(scaffold_num)),
    aes(x = (start + end) / 2 / 1e6, y = as.numeric(chr), label = scaffold_num),
    size = 3, angle = 0, vjust = 0.5, hjust = 0.5
  ) +
  scale_y_continuous(
    breaks = 1:length(levels(plot_data$chr)),
    labels = levels(plot_data$chr),
    name = "Chromosome"
  ) +
  scale_x_continuous(name = "Position (Mbp)") +
  scale_fill_identity() +
  labs(title = "Final Assembly Structure by Chromosome") +
  theme_minimal() +
  theme(panel.grid = element_blank())


## FINAL AGP
# Load AGP file
agp <- read.table("Pchalceus_SW.fasta.agp", header = FALSE, sep = "\t")

# Name columns according to AGP format
colnames(agp) <- c("chr", "start", "end", "part", "type", "scaffold", "scaf_start", "scaf_end", "orientation")

# Separate components (W) and gaps (U)
components <- agp %>% filter(type == "W")
gaps <- agp %>% filter(type == "U")

# Extract scaffold number for labeling
components$scaffold_num <- str_replace(components$scaffold, "scaffold_", "")

# Assign colors for chromosomes
chrom_colors <- scales::hue_pal()(length(unique(components$chr)))
names(chrom_colors) <- unique(components$chr)

# Add fill color based on chromosome and orientation
components <- components %>%
  mutate(
    chr_color = chrom_colors[chr],
    fill = ifelse(orientation == "+", chr_color, scales::muted(chr_color))  # darker shade for reverse
  )

# Gaps are light gray
gaps$fill <- "lightgray"

# Combine for plotting
plot_data <- bind_rows(
  components %>% select(chr, start, end, fill, scaffold_num),
  gaps %>% mutate(scaffold_num = NA) %>% select(chr, start, end, fill, scaffold_num)
)

# Order chromosomes
plot_data$chr <- factor(plot_data$chr, levels = paste0("CHR", 1:12))

ggplot(plot_data) +
  geom_rect(aes(xmin = start / 1e6, xmax = end / 1e6,
                ymin = as.numeric(chr) - 0.4, ymax = as.numeric(chr) + 0.4,
                fill = fill), color = "black") +
  geom_text(
    data = filter(plot_data, !is.na(scaffold_num)),
    aes(x = (start + end) / 2 / 1e6, y = as.numeric(chr), label = scaffold_num),
    size = 3, angle = 0, vjust = 0.5, hjust = 0.5
  ) +
  scale_y_continuous(
    breaks = 1:length(levels(plot_data$chr)),
    labels = levels(plot_data$chr),
    name = "Chromosome"
  ) +
  scale_x_continuous(name = "Position (Mbp)") +
  scale_fill_identity() +
  labs(title = "Final Assembly Structure by Chromosome") +
  theme_minimal() +
  theme(panel.grid = element_blank())
