setwd("~/Downloads")

library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

stats_SP <- read.csv("stats_SP.csv", header = TRUE)

filtered_stats_SP <- stats_SP[, c(1, 2, 3, 4, 20, 22)]

colnames(filtered_stats_SP)[colnames(filtered_stats_SP) == "FstWC_S_HUE_T_S_COT_S_egg"] <- "FstWC"
colnames(filtered_stats_SP)[colnames(filtered_stats_SP) == "lseff_S_HUE_T_S_COT_S_egg"] <- "lseff"

filtered_stats_SP <- subset(filtered_stats_SP, lseff > 400)

# Convert to numeric and handle warnings
filtered_stats_SP$FstWC <- as.numeric(as.character(filtered_stats_SP$FstWC))

# Check for any NA values introduced during conversion and remove them if necessary
filtered_stats_SP <- filtered_stats_SP[!is.na(filtered_stats_SP$FstWC), ]

# Replace values in FstWC column that are smaller than zero with zero
filtered_stats_SP$FstWC[filtered_stats_SP$FstWC < 0] <- 0
colfuncR <- colorRampPalette(c("black", "red"))
filtered_stats_SP$col <- colfuncR(100)[as.integer((filtered_stats_SP$FstWC/0.8)*100)+1]

# Function to find the first and last positions crossing the FstWC threshold for each scaffold
find_peak_positions <- function(df, threshold = 0.05, span = 0.05) {
  # Group by scaffold
  df %>%
    group_by(scaffold) %>%
    # Apply operations for each scaffold
    group_modify(~ {
      # Smooth the FstWC values using LOESS
      loess_fit <- loess(FstWC ~ mid, data = .x, span = span)
      smoothed_fst <- predict(loess_fit)
      
      # Add smoothed FstWC values to the dataframe
      .x <- mutate(.x, smoothed_FstWC = smoothed_fst)
      
      # Find the first and last positions crossing the threshold
      first_peak <- filter(.x, smoothed_FstWC > threshold) %>% slice_head(n = 1)
      last_peak <- filter(.x, smoothed_FstWC > threshold) %>% slice_tail(n = 1)
      
      # Combine the first and last peaks
      bind_rows(first_peak, last_peak)
    }) %>%
    ungroup()
}

# Apply the function to the filtered_stats_SP dataframe
threshold <- 0.1
span <- 0.05
peak_positions <- find_peak_positions(filtered_stats_SP, threshold, span)
peak_positions <- as.data.frame(peak_positions)

# Optionally save to a CSV file
write.csv(peak_positions, "peak_positions.csv", row.names = FALSE)

# Function to plot LOESS fit for a specific scaffold
plot_loess_fit <- function(df, scaffold, threshold = 0.05, span = 0.05) {
  # Filter data for the selected scaffold
  scaffold_data <- df %>% filter(scaffold == !!scaffold)
  
  # Check if data is empty
  if (nrow(scaffold_data) == 0) {
    stop(paste("Error: No data found for scaffold", scaffold))
  }
  
  # Fit LOESS
  loess_fit <- loess(FstWC ~ mid, data = scaffold_data, span = span)
  scaffold_data$smoothed_FstWC <- predict(loess_fit)
  
  # Create the plot
  p <- ggplot(scaffold_data, aes(x = mid)) +
    geom_point(aes(y = FstWC, color = FstWC), alpha = 0.7, size = 1.5) +  # Color based on FstWC
    geom_line(aes(y = smoothed_FstWC), color = "blue", linewidth = 1.2) +       # LOESS line in blue
    geom_hline(yintercept = threshold, linetype = "dashed", color = "black") +
    scale_color_gradient(low = "black", high = "red", name = "FstWC") +   # Gradient for FstWC
    labs(
      title = scaffold,
      x = "Position (mid)",
      y = "FstWC"
    ) +
    theme_minimal()
  
  return(p)
}

# Plot for a single chromosome (CM008230.1_RagTag, chromosome 1)
plot_loess_fit(filtered_stats_SP, scaffold = "CM008230.1_RagTag", threshold = 0.05, span = 0.05)

# List of scaffold names
scaffolds <- c(
  "CM008230.1_RagTag",
  "CM008231.1_RagTag",
  "CM008233.1_RagTag",
  "CM008234.1_RagTag",
  "CM008235.1_RagTag",
  "CM008236.1_RagTag",
  "CM008237.1_RagTag",
  "CM008238.1_RagTag",
  "CM008239.1_RagTag",
  "CM008240.1_RagTag"
)

plot_loess_fit(filtered_stats_SP, scaffold = "CM008240.1_RagTag", threshold = 0.05, span = 0.05)

# Generate a list of plots for each scaffold
plots <- lapply(scaffolds, function(scaffold) {
  plot_loess_fit(filtered_stats_SP, scaffold = scaffold, threshold = 0.05, span = 0.05)
})

# Arrange the plots in a 3x4 grid (3 rows and 4 columns)
grid.arrange(grobs = plots, ncol = 4, nrow = 3)
