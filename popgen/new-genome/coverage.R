setwd("/scratch/leuven/357/vsc35707/popgen/coverage")

library(ggplot2)
library(ggfortify)
library(dplyr)
library(stringr)
library(gridExtra)

GC129388 <- read.delim("GC129388.raw.pcov", header=FALSE, comment.char="#")
GC129389 <- read.delim("GC129389.raw.pcov", header=FALSE, comment.char="#")
GC129394 <- read.delim("GC129394.raw.pcov", header=FALSE, comment.char="#")
GC129395 <- read.delim("GC129395.raw.pcov", header=FALSE, comment.char="#")

names(GC129388) <- c("region", "coverage", "length")
names(GC129389) <- c("region", "coverage", "length")
names(GC129394) <- c("region", "coverage", "length")
names(GC129395) <- c("region", "coverage", "length")

lengths <- read.table("P_chalceus_final_lengths.txt", h=T)

# Input data
lengths <- data.frame(
  chromosome = paste0("CHR", 1:12),
  length = c(
    53807401, 53932439, 43748134, 44639843, 65108617, 45581523,
    42708394, 35503798, 27201578, 15106519, 21073294, 28338870
  )
)

# Calculate cumulative starts and ends
lengths$chromStarts <- cumsum(c(1, head(lengths$length, -1)))
lengths$chromEnds <- lengths$chromStarts + lengths$length - 1
lengths$chromMid <- (lengths$chromStarts + lengths$chromEnds) / 2

# Add scaffold-like names (e.g., scaffold_01)
lengths$scaffold <- sprintf("scaffold_%02d", 1:nrow(lengths))

# Rename columns
chrom_coords <- data.frame(
  chromosome = 1:nrow(lengths),
  scaffold = lengths$scaffold,
  length = lengths$length,
  chromLengths = lengths$length,
  chromStarts = lengths$chromStarts,
  chromEnds = lengths$chromEnds,
  chromMid = lengths$chromMid
)

# View result
print(chrom_coords)


##GC129388
# Filter out rows without ">" symbol in the "region" column
GC129388 <- GC129388 %>%
  filter(grepl(">", region))

# Remove ">" from region column if present
GC129388$region <- gsub(">", "", GC129388$region)

# Filter rows that start with "CHR"
GC129388 <- GC129388 %>%
  filter(str_starts(region, "CHR")) %>%
  mutate(chromosome = as.numeric(str_extract(region, "\\d+")))

# Remove everything before the last underscore in the "region" column
GC129388 <- GC129388 %>%
  mutate(region = gsub(".*_", "", region))

# Convert "region" and "length" columns to numeric (if they aren't already)
GC129388 <- GC129388 %>%
  mutate(
    region = as.numeric(region),
    length = as.numeric(length)
  ) %>%
  # Create a new column "position" by multiplying "region" and "length" columns
  mutate(position = region * length)

# Sort the data frame by the "chromosome" column
GC129388 <- GC129388 %>%
  arrange(chromosome) %>%
  group_by(chromosome)

GC129388 <- as.data.frame(GC129388)

#Merge data with scaffold lengths
merged_GC129388 <- merge(GC129388, chrom_coords, by = "chromosome", all.x = TRUE)

# Remove the length.y column from the merged data frame
merged_GC129388 <- merged_GC129388[, !(names(merged_GC129388) %in% c("length.y"))]

# Rename the position column to match the original column name
names(merged_GC129388)[names(merged_GC129388) == "length.x"] <- "window_length"

#Change names
merged_GC129388$genomic_position <- (merged_GC129388$position + merged_GC129388$chromStarts)-1

#Add new columns
merged_GC129388$individual <- "GC129388"
merged_GC129388$population <- "Nieuwpoort"

average_coverage <- mean(merged_GC129388$coverage)
merged_GC129388 <- transform(merged_GC129388, normalized_coverage = coverage / average_coverage)

ggplot(merged_GC129388, aes(x=genomic_position, y = normalized_coverage))+ geom_point()+ theme_minimal()


# Define sample-to-population mapping
sample_info <- data.frame(
  sample = c("GC136107", "GC136108", "GC136109", "GC136110", "GC136111",
             "GC136112", "GC136113", "GC136114", "GC136115", "GC136116"),
  population = c("S_HUE_T", "S_HUE_T", "S_COT_S", "S_COT_S", "S_COT_S",
                 "S_COT_S", "S_COT_S", "S_HUE_T", "S_HUE_T", "S_HUE_T"),
  stringsAsFactors = FALSE
)

# Create list to store all processed data frames
processed_list <- list()

# Loop over each sample
for (i in 1:nrow(sample_info)) {
  sample_id <- sample_info$sample[i]
  pop <- sample_info$population[i]
  
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
      chromosome = as.numeric(str_extract(region, "\\d+")),
      region = gsub(".*_", "", region),
      region = as.numeric(region),
      length = as.numeric(length),
      position = region * length
    ) %>%
    arrange(chromosome) %>%
    group_by(chromosome) %>%
    ungroup()
  
  # Merge with chromosome coordinates
  merged_df <- merge(data, chrom_coords, by = "chromosome", all.x = TRUE)
  merged_df <- merged_df[, !(names(merged_df) %in% c("length.y"))]
  names(merged_df)[names(merged_df) == "length.x"] <- "window_length"
  
  # Compute genomic position and add metadata
  merged_df$genomic_position <- merged_df$position + merged_df$chromStarts - 1
  merged_df$individual <- sample_id
  merged_df$population <- pop
  average_coverage <- mean(merged_df$coverage, na.rm = TRUE)
  merged_df$normalized_coverage <- merged_df$coverage / average_coverage
  
  # Save each dataframe as a separate variable in the global environment
  assign(sample_id, merged_df, envir = .GlobalEnv)
  
  # Also save in a list if needed
  processed_list[[sample_id]] <- merged_df
}


merged_data <- do.call(rbind, processed_list)

ggplot(merged_data, aes(x = genomic_position, y = normalized_coverage, group = individual, color = population)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~ chromosome, scales = "free_x", ncol = 3) +
  scale_color_manual(values = c("S_HUE_T" = "blue", "S_COT_S" = "red")) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    legend.title = element_blank()
  ) +
  labs(
    x = "Genomic position (bp)",
    y = "Normalized coverage",
    title = "Normalized coverage per chromosome by individual and population"
  )

# Create an empty list to hold the plots
plots <- list()

# Get all chromosomes present in the data
chromosomes <- sort(unique(merged_data$chromosome))

# Generate plots
for (chr in chromosomes) {
  chr_data <- merged_data %>% filter(chromosome == chr)
  
  p <- ggplot(chr_data, aes(x = genomic_position / 1e6,  # Mbp
                            y = normalized_coverage,
                            group = individual,
                            color = population)) +
    geom_line(alpha = 0.3) +  # More transparent lines
    scale_color_manual(values = c("S_HUE_T" = "blue", "S_COT_S" = "red")) +
    theme_minimal() +
    labs(
      title = paste("Chromosome", chr),
      x = "Genomic position (Mbp)",
      y = "Normalized coverage"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 10, face = "bold")
    )
  
  plots[[length(plots) + 1]] <- p
}
grid.arrange(grobs = plots, ncol = 4)

# Define sample-to-sex mapping
sex_info <- data.frame(
  individual = c("GC136107", "GC136108", "GC136109", "GC136110", "GC136111",
                 "GC136112", "GC136113", "GC136114", "GC136115", "GC136116"),
  sex = c("F", "F", "F", "M", "F", "F", "F", "F", "F", "M"),
  stringsAsFactors = FALSE
)
merged_data <- merged_data %>%
  left_join(sex_info, by = "individual")

# Filter for chromosomes 10, 11 and 12
sex_chrom_data <- merged_data %>%
  filter(chromosome %in% c(10, 11, 12))

# Define sex color mapping
sex_colors <- c("F" = "pink", "M" = "lightblue")

# Generate plot list
plots_sex <- list()

for (chr in c(10, 11, 12)) {
  chr_data <- sex_chrom_data %>% filter(chromosome == chr)
  
  p <- ggplot(chr_data, aes(x = genomic_position / 1e6,
                            y = normalized_coverage,
                            group = individual,
                            color = sex)) +
    geom_line(alpha = 0.8, size = 0.8) +
    scale_color_manual(values = sex_colors) +
    theme_minimal() +
    labs(
      title = paste("Chromosome", chr),
      x = "Genomic position (Mbp)",
      y = "Normalized coverage"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12, face = "bold")
    )
  
  plots_sex[[length(plots_sex) + 1]] <- p
}

# Display the two plots side by side
grid.arrange(grobs = plots_sex, ncol = 3)
