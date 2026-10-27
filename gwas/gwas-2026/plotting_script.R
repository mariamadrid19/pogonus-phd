
# =========================================================
# GWAS on the new LW genome
# Manhattan plots — Nieuwpoort and Dudzele 2026
# =========================================================

library(data.table)
library(dplyr)
library(ggplot2)

# =========================================================
# 1. Settings
# =========================================================

setwd("~/Documents/phd-pogonus/GWAS-2026-02")

# =========================================================
# 2. Calculate chromosome coordinates
# =========================================================

chrom.coords <- function(scafL, chromNames, gap = 9000000) {
  
  chromosome <- vector()
  chromLengths <- vector()
  chromStarts <- vector()
  chromEnds <- vector()
  chromMid <- vector()
  
  chrom <- 1
  endLast <- 0
  
  scafCurrent <- subset(
    scafL,
    chromosome == chromNames[1]
  )
  
  chromosome[chrom] <- chrom
  chromLengths[chrom] <- sum(scafCurrent$length)
  chromStarts[chrom] <- endLast + 1
  chromEnds[chrom] <- endLast + chromLengths[chrom]
  chromMid[chrom] <- endLast + chromLengths[chrom] / 2
  
  endLast <- chromEnds[chrom]
  chrom <- chrom + 1
  
  for (i in 2:length(chromNames)) {
    
    chromosome[chrom] <- chrom
    
    scafCurrent <- subset(
      scafL,
      chromosome == chromNames[i]
    )
    
    chromLengths[chrom] <- sum(scafCurrent$length)
    chromStarts[chrom] <- endLast + gap + 1
    chromEnds[chrom] <- endLast + gap + chromLengths[chrom]
    chromMid[chrom] <- endLast + gap + chromLengths[chrom] / 2
    
    endLast <- chromEnds[chrom]
    chrom <- chrom + 1
  }
  
  table <- as.data.frame(
    cbind(
      chromosome,
      chromLengths,
      chromStarts,
      chromEnds,
      chromMid
    )
  )
  
  return(table)
}

# ---------------------------------------------------------
# Read chromosome lengths
# ---------------------------------------------------------

scafL <- read.table(
  "Pchalceus_LW_final.fasta.sorted.lengths.txt",
  header = FALSE,
  stringsAsFactors = FALSE
)

colnames(scafL) <- c("chrom_name", "length")

# Extract CHR entries only
scafL <- scafL[
  grepl("^CHR", scafL$chrom_name),
]

# Extract chromosome numbers
scafL$chrom_id <- as.numeric(
  gsub("CHR", "", scafL$chrom_name)
)

scafL$chromosome <- scafL$chrom_id

scafL$length <- as.numeric(scafL$length)

# Calculate coordinates for 11 chromosomes
chrom_coords <- chrom.coords(
  scafL = scafL,
  chromNames = 1:11
)

print(chrom_coords)

# Merge chromosome coordinates
scaf_coords <- merge(
  scafL,
  chrom_coords,
  by = "chromosome",
  all.x = TRUE
)

names(scaf_coords)[
  names(scaf_coords) == "chrom_name"
] <- "scaffold"

scaf_coords$scaffold <- as.character(
  scaf_coords$scaffold
)

scaf_coords <- scaf_coords %>%
  mutate(chromosome = as.character(chromosome))

# =========================================================
# 3. Datasets
# =========================================================

datasets <- list(
  
  Nieuwpoort = list(
    
    results_dir = "gemma-results-Nieuwpoort-LW-genome",
    
    output_dir = "manhattan-plots-Nieuwpoort",
    
    n_individuals = 192,
    
    phenotypes = c(
      "WL",
      "EL",
      "relMRWS",
      "time5trials",
      "time3trials",
      "emerged95s5trials",
      "emerged95s3trials",
      "endPosition3trials",
      "action5trials",
      "action3trials"
    ),
    
    phenotype_titles = c(
      WL = "Wing length",
      EL = "Elytra length",
      relMRWS = "Relative MRWS",
      time5trials = "Time (5 trials)",
      time3trials = "Time (3 trials)",
      emerged95s5trials = "Emerged 95s (5 trials)",
      emerged95s3trials = "Emerged 95s (3 trials)",
      endPosition3trials = "End position (3 trials)",
      action5trials = "Action (5 trials)",
      action3trials = "Action (3 trials)"
    )
  ),
  
  Dudzele = list(
    
    results_dir = "gemma-results-Dudzele-LW-genome",
    
    output_dir = "manhattan-plots-Dudzele",
    
    n_individuals = 180,
    
    phenotypes = c(
      "WL",
      "EL",
      "time3trials",
      "endPosition3trials_binary",
      "relMRWS"
    ),
    
    phenotype_titles = c(
      WL = "Wing length",
      EL = "Elytra length",
      time3trials = "Time (3 trials)",
      endPosition3trials_binary = "End position (3 trials)",
      relMRWS = "Relative MRWS"
    )
  )
)

# =========================================================
# 4. Loop over datasets
# =========================================================

for (POP in names(datasets)) {
  
  settings <- datasets[[POP]]
  
  results_dir <- settings$results_dir
  output_dir <- settings$output_dir
  phenotypes <- settings$phenotypes
  phenotype_titles <- settings$phenotype_titles
  n_individuals <- settings$n_individuals
  
  dir.create(
    output_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )
  
  cat("\n====================================\n")
  cat("Population:", POP, "\n")
  cat("====================================\n")
  
  # =======================================================
  # 5. Loop over all phenotypes
  # =======================================================
  
  for (PHENO in phenotypes) {
    
    cat("\n====================================\n")
    cat("Processing:", PHENO, "\n")
    cat("====================================\n")
    
    # -----------------------------------------------------
    # Read GEMMA results
    # -----------------------------------------------------
    
    file_path <- file.path(
      results_dir,
      PHENO,
      paste0("gemma_lmm_", PHENO, ".assoc.txt")
    )
    
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path)
      next
    }
    
    gemma_results <- fread(file_path)
    
    gemma_results <- gemma_results %>%
      mutate(
        chr = as.character(chr),
        p_wald = as.numeric(p_wald),
        ps = as.numeric(ps)
      ) %>%
      left_join(
        scaf_coords,
        by = c("chr" = "chromosome")
      ) %>%
      mutate(
        pos_cumulative = ps + chromStarts
      )
    
    # Remove invalid positions and p-values
    gemma_results <- gemma_results %>%
      filter(
        !is.na(pos_cumulative),
        !is.na(p_wald),
        p_wald > 0,
        p_wald <= 1
      )
    
    if (nrow(gemma_results) == 0) {
      warning("No valid SNPs for ", POP, " ", PHENO)
      next
    }
    
    # -----------------------------------------------------
    # Bonferroni threshold
    # -----------------------------------------------------
    
    num_tests <- nrow(gemma_results)
    
    bonferroni_threshold <- 0.05 / num_tests
    
    bonferroni_threshold_log <- -log10(
      bonferroni_threshold
    )
    
    # -----------------------------------------------------
    # Dynamic scale bar position
    # -----------------------------------------------------
    
    # Determine height from the actual GWAS results.
    # Ensure that both significance thresholds are included.
    
    plot_height <- max(
      -log10(gemma_results$p_wald),
      bonferroni_threshold_log + 1,
      na.rm = TRUE
    )
    
    # Position scale bar relative to plot height
    scale_bar_y <- plot_height * 0.12
    scale_label_y <- plot_height * 0.07
    
    # -----------------------------------------------------
    # Chromosome colors
    # -----------------------------------------------------
    
    chromosome_order <- scaf_coords$chromosome
    
    gemma_results$chr <- factor(
      gemma_results$chr,
      levels = chromosome_order
    )
    
    chrom_colors <- setNames(
      rep(
        c("lightgray", "darkgray"),
        length.out = length(chromosome_order)
      ),
      chromosome_order
    )
    
    # -----------------------------------------------------
    # Number of individuals and plot title
    # -----------------------------------------------------
    
    plot_title <- paste0(
      phenotype_titles[[PHENO]],
      " (n=", n_individuals, ") ",
      POP,
      " 2026"
    )
    
    # -----------------------------------------------------
    # Manhattan plot
    # -----------------------------------------------------
    
    manhattan_pwald <- ggplot(
      gemma_results,
      aes(
        x = pos_cumulative,
        y = -log10(p_wald)
      )
    ) +
      
      geom_point(
        data = subset(
          gemma_results,
          p_wald >= bonferroni_threshold
        ),
        aes(color = chr),
        size = 2,
        alpha = 1
      ) +
      
      geom_point(
        data = subset(
          gemma_results,
          p_wald < bonferroni_threshold
        ),
        color = "red",
        size = 2,
        alpha = 1
      ) +
      
      geom_hline(
        yintercept = bonferroni_threshold_log,
        color = "red",
        linetype = "dashed",
        linewidth = 0.5
      ) +
      
      geom_hline(
        yintercept = -log10(1e-5),
        color = "blue",
        linetype = "dashed",
        linewidth = 0.5
      ) +
      
      annotate(
        "text",
        x = max(gemma_results$pos_cumulative) * 0.88,
        y = bonferroni_threshold_log + 0.5,
        label = "Bonferroni threshold",
        color = "red",
        size = 4
      ) +
      
      annotate(
        "text",
        x = max(gemma_results$pos_cumulative) * 0.873,
        y = -log10(1e-5) - 0.5,
        label = "Suggestive threshold",
        color = "blue",
        size = 4
      ) +
      
      scale_color_manual(
        values = chrom_colors
      ) +
      
      scale_x_continuous(
        breaks = scaf_coords$chromMid,
        labels = scaf_coords$chromosome
      ) +
      
      # ---------------------------------------------------
    # Dynamic 20 Mb scale bar
    # ---------------------------------------------------
    
    annotate(
      "segment",
      x = 5e6,
      xend = 25e6,
      y = scale_bar_y,
      yend = scale_bar_y,
      linewidth = 0.5
    ) +
      
      annotate(
        "text",
        x = 15e6,
        y = scale_label_y,
        label = "20 Mb",
        size = 3.5
      ) +
      
      theme_minimal() +
      
      labs(
        title = plot_title,
        subtitle = "Sex included as covariate",
        x = NULL,
        y = expression(-log[10](p))
      ) +
      
      theme(
        axis.text.x = element_text(hjust = 0.5),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(
          size = 16,
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          size = 12,
          hjust = 0.5
        ),
        axis.title = element_text(size = 16)
      )
    
    # -----------------------------------------------------
    # Save plot
    # -----------------------------------------------------
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0("manhattan_", PHENO, ".pdf")
      ),
      plot = manhattan_pwald,
      width = 14,
      height = 6
    )
    
    ggsave(
      filename = file.path(
        output_dir,
        paste0("manhattan_", PHENO, ".png")
      ),
      plot = manhattan_pwald,
      width = 14,
      height = 6,
      dpi = 300
    )
    
    # -----------------------------------------------------
    # Summary
    # -----------------------------------------------------
    
    cat("Number of individuals:", n_individuals, "\n")
    
    cat("Number of SNPs:", num_tests, "\n")
    
    cat(
      "Bonferroni threshold:",
      bonferroni_threshold,
      "\n"
    )
    
    cat(
      "Significant SNPs:",
      sum(
        gemma_results$p_wald < bonferroni_threshold
      ),
      "\n"
    )
    
    cat(
      "Saved Manhattan plot for:",
      POP,
      PHENO,
      "\n"
    )
    
  }
  
  cat("\nCompleted Manhattan plots for:", POP, "\n")
}

cat("\nAll Nieuwpoort and Dudzele Manhattan plots completed!\n")



# =========================================================
# 6. Combined Belgium GWAS — Manhattan plots
# =========================================================

results_dir <- "combined/gemma-results"
output_dir <- "combined/manhattan-plots"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

n_individuals <- 372

phenotypes <- c(
  "WL",
  "EL",
  "time3trials",
  "relMRWS"
)

phenotype_titles <- c(
  WL = "Wing length",
  EL = "Elytra length",
  time3trials = "Time (3 trials)",
  relMRWS = "Relative MRWS"
)

# =========================================================
# 7. Loop over all combined phenotypes
# =========================================================

for (PHENO in phenotypes) {
  
  cat("\n====================================\n")
  cat("Processing combined Belgium:", PHENO, "\n")
  cat("====================================\n")
  
  # -------------------------------------------------------
  # Read GEMMA results
  # -------------------------------------------------------
  
  file_path <- file.path(
    results_dir,
    PHENO,
    paste0("gemma_lmm_", PHENO, ".assoc.txt")
  )
  
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    next
  }
  
  gemma_results <- fread(file_path)
  
  gemma_results <- gemma_results %>%
    mutate(
      chr = as.character(chr),
      p_wald = as.numeric(p_wald),
      ps = as.numeric(ps)
    ) %>%
    left_join(
      scaf_coords,
      by = c("chr" = "chromosome")
    ) %>%
    mutate(
      pos_cumulative = ps + chromStarts
    )
  
  # Remove invalid positions and p-values
  gemma_results <- gemma_results %>%
    filter(
      !is.na(pos_cumulative),
      !is.na(p_wald),
      p_wald > 0,
      p_wald <= 1
    )
  
  if (nrow(gemma_results) == 0) {
    warning("No valid SNPs for combined Belgium ", PHENO)
    next
  }
  
  # -------------------------------------------------------
  # Bonferroni threshold
  # -------------------------------------------------------
  
  num_tests <- nrow(gemma_results)
  
  bonferroni_threshold <- 0.05 / num_tests
  
  bonferroni_threshold_log <- -log10(
    bonferroni_threshold
  )
  
  # -------------------------------------------------------
  # Dynamic scale bar position
  # -------------------------------------------------------
  
  plot_height <- max(
    -log10(gemma_results$p_wald),
    bonferroni_threshold_log + 1,
    na.rm = TRUE
  )
  
  scale_bar_y <- plot_height * 0.12
  scale_label_y <- plot_height * 0.07
  
  # -------------------------------------------------------
  # Chromosome colors
  # -------------------------------------------------------
  
  chromosome_order <- scaf_coords$chromosome
  
  gemma_results$chr <- factor(
    gemma_results$chr,
    levels = chromosome_order
  )
  
  chrom_colors <- setNames(
    rep(
      c("lightgray", "darkgray"),
      length.out = length(chromosome_order)
    ),
    chromosome_order
  )
  
  # -------------------------------------------------------
  # Plot title
  # -------------------------------------------------------
  
  plot_title <- paste0(
    phenotype_titles[[PHENO]],
    " (n=", n_individuals, ") Belgium 2026"
  )
  
  # -------------------------------------------------------
  # Manhattan plot
  # -------------------------------------------------------
  
  manhattan_pwald <- ggplot(
    gemma_results,
    aes(
      x = pos_cumulative,
      y = -log10(p_wald)
    )
  ) +
    
    geom_point(
      data = subset(
        gemma_results,
        p_wald >= bonferroni_threshold
      ),
      aes(color = chr),
      size = 2,
      alpha = 1
    ) +
    
    geom_point(
      data = subset(
        gemma_results,
        p_wald < bonferroni_threshold
      ),
      color = "red",
      size = 2,
      alpha = 1
    ) +
    
    geom_hline(
      yintercept = bonferroni_threshold_log,
      color = "red",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    geom_hline(
      yintercept = -log10(1e-5),
      color = "blue",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    
    annotate(
      "text",
      x = max(gemma_results$pos_cumulative) * 0.88,
      y = bonferroni_threshold_log + 0.5,
      label = "Bonferroni threshold",
      color = "red",
      size = 4
    ) +
    
    annotate(
      "text",
      x = max(gemma_results$pos_cumulative) * 0.873,
      y = -log10(1e-5) - 0.5,
      label = "Suggestive threshold",
      color = "blue",
      size = 4
    ) +
    
    scale_color_manual(
      values = chrom_colors
    ) +
    
    scale_x_continuous(
      breaks = scaf_coords$chromMid,
      labels = scaf_coords$chromosome
    ) +
    
    # -----------------------------------------------------
  # Dynamic 20 Mb scale bar
  # -----------------------------------------------------
  
  annotate(
    "segment",
    x = 5e6,
    xend = 25e6,
    y = scale_bar_y,
    yend = scale_bar_y,
    linewidth = 0.5
  ) +
    
    annotate(
      "text",
      x = 15e6,
      y = scale_label_y,
      label = "20 Mb",
      size = 3.5
    ) +
    
    theme_minimal() +
    
    labs(
      title = plot_title,
      subtitle = "Sex and population included as covariates",
      x = NULL,
      y = expression(-log[10](p))
    ) +
    
    theme(
      axis.text.x = element_text(hjust = 0.5),
      axis.ticks.x = element_blank(),
      legend.position = "none",
      plot.title = element_text(
        size = 16,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 12,
        hjust = 0.5
      ),
      axis.title = element_text(size = 16)
    )
  
  # -------------------------------------------------------
  # Save plot
  # -------------------------------------------------------
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0("manhattan_", PHENO, ".pdf")
    ),
    plot = manhattan_pwald,
    width = 14,
    height = 6
  )
  
  ggsave(
    filename = file.path(
      output_dir,
      paste0("manhattan_", PHENO, ".png")
    ),
    plot = manhattan_pwald,
    width = 14,
    height = 6,
    dpi = 300
  )
  
  # -------------------------------------------------------
  # Summary
  # -------------------------------------------------------
  
  cat("Number of individuals:", n_individuals, "\n")
  cat("Number of SNPs:", num_tests, "\n")
  cat("Bonferroni threshold:", bonferroni_threshold, "\n")
  
  cat(
    "Significant SNPs:",
    sum(gemma_results$p_wald < bonferroni_threshold),
    "\n"
  )
  
  cat(
    "Saved combined Belgium Manhattan plot for:",
    PHENO,
    "\n"
  )
}

cat("\nAll combined Belgium Manhattan plots completed!\n")
