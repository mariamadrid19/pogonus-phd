library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(purrr)

setwd("filtering-directory")

prefix <- "map12"
output_prefix <- paste0(prefix, ".clean")
n_chr <- 11
top_n_scaffolds <- 30

marker_files <- sprintf("%s_chr%d_m.clean.input", prefix, 1:n_chr)
agp_files <- sprintf("clean_agp/%s_chr%d.clean.agp", prefix, 1:n_chr)

read_marker_file <- function(f) {
  chr <- str_extract(f, "(?<=chr)\\d+")
  
  read_table(
    f,
    col_names = c("scaffold", "pos", "LG", "male_cM", "female_cM"),
    col_types = cols(
      scaffold = col_character(),
      pos = col_double(),
      LG = col_integer(),
      male_cM = col_double(),
      female_cM = col_double()
    )
  ) %>%
    mutate(chromosome = as.integer(chr)) %>%
    arrange(female_cM, scaffold, pos) %>%
    mutate(marker_rank = row_number())
}

markers <- map_dfr(marker_files, read_marker_file)

read_agp_file <- function(f) {
  chr <- str_extract(f, "(?<=chr)\\d+")
  
  df <- read_tsv(
    f,
    comment = "#",
    col_names = FALSE,
    col_types = cols(.default = col_character()),
    show_col_types = FALSE
  )
  
  colnames(df)[1:9] <- c(
    "object", "object_beg", "object_end", "part_number", "component_type",
    "component_id", "component_beg", "component_end", "orientation"
  )
  
  df %>%
    mutate(
      chromosome = as.integer(chr),
      chromosome_num = chromosome,
      object_beg = as.numeric(object_beg),
      object_end = as.numeric(object_end)
    )
}

agp <- map_dfr(agp_files, read_agp_file)

# p3
agp_scaff_summary <- agp %>%
  filter(component_type == "W") %>%
  group_by(chromosome, component_id) %>%
  summarise(
    start = min(object_beg),
    end = max(object_end),
    orientation = first(orientation),
    .groups = "drop"
  ) %>%
  arrange(chromosome, start) %>%
  group_by(chromosome) %>%
  mutate(scaffold_rank = row_number()) %>%
  ungroup()

p3 <- ggplot(
  agp_scaff_summary,
  aes(x = start, xend = end, y = scaffold_rank, yend = scaffold_rank, color = orientation)
) +
  geom_segment(linewidth = 2) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 3) +
  labs(
    title = "Cleaned anchored scaffold layout from AGP",
    x = "Anchored chromosome position (bp)",
    y = "Scaffold rank along chromosome",
    color = "Orientation"
  ) +
  theme_bw()

ggsave(paste0(output_prefix, "_p3_clean_AGP_layout.pdf"), p3, width = 14, height = 10)
ggsave(paste0(output_prefix, "_p3_clean_AGP_layout.png"), p3, width = 14, height = 10, dpi = 300)

# p5
markers_plot <- markers %>%
  group_by(chromosome, scaffold) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(chromosome) %>%
  slice_max(order_by = n, n = top_n_scaffolds, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(keep = TRUE) %>%
  right_join(markers, by = c("chromosome", "scaffold")) %>%
  mutate(scaffold_plot = ifelse(is.na(keep), "other", scaffold)) %>%
  select(-n, -keep)

p5 <- ggplot(markers_plot, aes(x = female_cM, y = pos, color = scaffold_plot)) +
  geom_point(size = 0.5, alpha = 0.8) +
  facet_wrap(~ chromosome, scales = "free_y", ncol = 3) +
  labs(
    title = "Cleaned physical marker positions vs genetic positions",
    x = "Female map position (cM)",
    y = "Marker position on scaffold (bp)",
    color = "Scaffold"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

ggsave(paste0(output_prefix, "_p5_clean_marker_positions.pdf"), p5, width = 14, height = 10)
ggsave(paste0(output_prefix, "_p5_clean_marker_positions.png"), p5, width = 14, height = 10, dpi = 300)

# p_agp
agp_scaff_blocks <- agp %>%
  filter(component_type == "W") %>%
  mutate(
    chromosome = paste0("CHR", chromosome_num)
  ) %>%
  arrange(chromosome_num, object_beg)

scaffold_labels <- agp_scaff_blocks %>%
  distinct(component_id) %>%
  mutate(scaffold_label = readr::parse_number(component_id))

agp_scaff_blocks <- agp_scaff_blocks %>%
  left_join(scaffold_labels, by = "component_id")

chr_levels <- paste0("CHR", n_chr:1)
agp_scaff_blocks$chromosome <- factor(agp_scaff_blocks$chromosome, levels = chr_levels)

p_agp <- ggplot(agp_scaff_blocks) +
  geom_rect(
    aes(
      xmin = object_beg / 1e6,
      xmax = object_end / 1e6,
      ymin = as.numeric(chromosome) - 0.4,
      ymax = as.numeric(chromosome) + 0.4,
      fill = factor(scaffold_label)
    ),
    color = "black",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      x = (object_beg + object_end) / 2e6,
      y = as.numeric(chromosome),
      label = scaffold_label
    ),
    size = 4
  ) +
  scale_y_continuous(
    breaks = seq_along(chr_levels),
    labels = chr_levels
  ) +
  labs(
    title = "Cleaned Final Assembly Structure by Chromosome",
    x = "Position (Mbp)",
    y = "Chromosome"
  ) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave(paste0(output_prefix, "_pAGP_clean_final_structure.pdf"), p_agp, width = 11, height = 11)
ggsave(paste0(output_prefix, "_pAGP_clean_final_structure.png"), p_agp, width = 11, height = 11, dpi = 300)

cat("Done. Plotted newest cleaned run only.\n")
