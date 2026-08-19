library(SVbyEye)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

# ============================================================
# Input paths
# ============================================================

paf.file <- "/Users/mariamadrid/Downloads/PchalceusLWrefV2.paf.gz"
fst.directory <- "/Users/mariamadrid/Downloads/LW"

chromosome.length.file <- file.path(
  fst.directory,
  "scaflengthsLW.txt"
)

chr.order <- paste0("CHR", 1:11)

# ============================================================
# 1. Import the PAF alignments
# Creates: paf.table
# ============================================================

if (!file.exists(paf.file)) {
  stop("PAF file does not exist: ", paf.file)
}

paf.table <- SVbyEye::readPaf(
  paf.file = paf.file,
  include.paf.tags = FALSE
)

if (nrow(paf.table) == 0) {
  stop("The imported PAF table contains no alignments.")
}

required.paf.columns <- c(
  "q.name",
  "q.len",
  "q.start",
  "q.end",
  "strand",
  "t.name",
  "t.len",
  "t.start",
  "t.end",
  "n.match"
)

missing.paf.columns <- setdiff(
  required.paf.columns,
  names(paf.table)
)

if (length(missing.paf.columns) > 0) {
  stop(
    "The following required PAF columns are missing: ",
    paste(missing.paf.columns, collapse = ", ")
  )
}

message(
  "Imported ",
  format(nrow(paf.table), big.mark = ","),
  " PAF alignments."
)

# Confirm that LW chromosomes are present as PAF targets
missing.target.chromosomes <- setdiff(
  chr.order,
  unique(as.character(paf.table$t.name))
)

if (length(missing.target.chromosomes) > 0) {
  warning(
    "The following chromosomes were not found in paf.table$t.name: ",
    paste(missing.target.chromosomes, collapse = ", "),
    ". Confirm that LW is the minimap2 target/reference assembly."
  )
}

# ============================================================
# 2. Read LW chromosome lengths
# Creates: chrom_data_chr
# ============================================================

if (!file.exists(chromosome.length.file)) {
  stop(
    "Chromosome-length file does not exist: ",
    chromosome.length.file
  )
}

chrom_data <- read.table(
  chromosome.length.file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required.length.columns <- c("chrom", "length")

missing.length.columns <- setdiff(
  required.length.columns,
  names(chrom_data)
)

if (length(missing.length.columns) > 0) {
  stop(
    "The chromosome-length file is missing: ",
    paste(missing.length.columns, collapse = ", "),
    ". Available columns are: ",
    paste(names(chrom_data), collapse = ", ")
  )
}

chrom_data_chr <- chrom_data %>%
  transmute(
    chrom = as.character(chrom),
    length = suppressWarnings(
      as.numeric(as.character(length))
    )
  ) %>%
  filter(chrom %in% chr.order) %>%
  mutate(
    chromosome = as.integer(
      sub("^CHR", "", chrom)
    )
  ) %>%
  arrange(chromosome)

if (nrow(chrom_data_chr) != 11) {
  stop(
    "Expected 11 main chromosomes in ",
    chromosome.length.file,
    ", but found ",
    nrow(chrom_data_chr),
    ". Chromosomes found: ",
    paste(chrom_data_chr$chrom, collapse = ", ")
  )
}

if (
  anyNA(chrom_data_chr$length) ||
  any(chrom_data_chr$length <= 0)
) {
  stop(
    "Some chromosome lengths are missing, nonnumeric, or nonpositive."
  )
}

message("Imported lengths for all 11 LW chromosomes.")

# ============================================================
# 3. Identify the best SW chromosome for each LW chromosome
# Creates: best.pairs
# ============================================================

pair.summary <- paf.table %>%
  mutate(
    t.name = as.character(t.name),
    q.name = as.character(q.name),
    n.match = suppressWarnings(
      as.numeric(as.character(n.match))
    )
  ) %>%
  filter(
    t.name %in% chr.order,
    q.name %in% chr.order,
    !is.na(n.match)
  ) %>%
  group_by(t.name, q.name) %>%
  summarise(
    matched.bp = sum(n.match, na.rm = TRUE),
    alignment.count = n(),
    .groups = "drop"
  )

if (nrow(pair.summary) == 0) {
  stop(
    "No chromosome-to-chromosome alignments were found. ",
    "Check the chromosome names in q.name and t.name."
  )
}

best.pairs <- pair.summary %>%
  group_by(t.name) %>%
  slice_max(
    order_by = matched.bp,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    chromosome.number = as.integer(
      sub("^CHR", "", t.name)
    )
  ) %>%
  arrange(chromosome.number) %>%
  select(
    t.name,
    q.name,
    matched.bp,
    alignment.count
  )

missing.best.pairs <- setdiff(
  chr.order,
  best.pairs$t.name
)

if (length(missing.best.pairs) > 0) {
  warning(
    "No best SW chromosome match was found for: ",
    paste(missing.best.pairs, collapse = ", ")
  )
}

message("Best LW–SW chromosome matches:")

print(best.pairs)

# ============================================================
# 4. Locate the four sets of FST files
# ============================================================

if (!dir.exists(fst.directory)) {
  stop(
    "FST directory does not exist: ",
    fst.directory
  )
}

# The order here determines the order of the FST tracks
population.info <- data.frame(
  population = c(
    "Belgium",
    "France",
    "Portugal",
    "Spain"
  ),
  population.number = 0:3,
  comparison = c(
    "Belgium SW/LW",
    "France SW/LW",
    "Portugal SW/LW",
    "Spain SW/LW"
  ),
  stringsAsFactors = FALSE
)

# Creates: comparison_labels
comparison_labels <- population.info$comparison

sort_stats_files <- function(files) {
  
  chromosome.number <- as.integer(
    sub(
      ".*\\.chr_([0-9]+)\\.pop[0-9]+\\.stats$",
      "\\1",
      basename(files)
    )
  )
  
  files[order(chromosome.number)]
}

file_lists <- vector(
  mode = "list",
  length = nrow(population.info)
)

names(file_lists) <- population.info$population

for (i in seq_len(nrow(population.info))) {
  
  pop.number <- population.info$population.number[i]
  
  file.pattern <- paste0(
    "^Pchal_LW\\.chr_",
    "(1|2|3|4|5|6|7|8|9|10|11)",
    "\\.pop",
    pop.number,
    "\\.stats$"
  )
  
  files <- list.files(
    path = fst.directory,
    pattern = file.pattern,
    full.names = TRUE
  )
  
  files <- sort_stats_files(files)
  
  if (length(files) != 11) {
    stop(
      "Expected 11 FST files for ",
      population.info$population[i],
      " (pop",
      pop.number,
      "), but found ",
      length(files),
      ".\nFiles found:\n",
      paste(files, collapse = "\n")
    )
  }
  
  # Check that each chromosome occurs exactly once
  detected.chromosomes <- as.integer(
    sub(
      ".*\\.chr_([0-9]+)\\.pop[0-9]+\\.stats$",
      "\\1",
      basename(files)
    )
  )
  
  if (!setequal(detected.chromosomes, 1:11)) {
    stop(
      "The files for ",
      population.info$population[i],
      " do not contain exactly chromosomes 1-11. ",
      "Detected: ",
      paste(detected.chromosomes, collapse = ", ")
    )
  }
  
  file_lists[[i]] <- files
  
  message(
    population.info$population[i],
    ": found all 11 FST files."
  )
}

# ============================================================
# 5. Import the FST files
# Creates: comp
# ============================================================

read_stats_files <- function(files) {
  
  imported.tables <- lapply(
    files,
    function(current.file) {
      
      x <- data.table::fread(
        current.file,
        data.table = FALSE
      )
      
      if (nrow(x) == 0) {
        warning(
          "The following statistics file is empty: ",
          current.file
        )
      }
      
      x$source_file <- basename(current.file)
      
      x
    }
  )
  
  result <- data.table::rbindlist(
    imported.tables,
    use.names = TRUE,
    fill = TRUE
  )
  
  result <- as.data.frame(result)
  
  if (nrow(result) == 0) {
    stop("The imported FST table contains no rows.")
  }
  
  if (!"scaffold" %in% names(result)) {
    stop(
      "No 'scaffold' column was found in the FST files. ",
      "Available columns are: ",
      paste(names(result), collapse = ", ")
    )
  }
  
  result$scaffold <- as.character(result$scaffold)
  
  result
}

comp <- lapply(
  file_lists,
  read_stats_files
)

names(comp) <- comparison_labels

for (i in seq_along(comp)) {
  message(
    comparison_labels[i],
    ": imported ",
    format(nrow(comp[[i]]), big.mark = ","),
    " FST windows."
  )
}

# ============================================================
# 6. Helper for finding FST and site-count columns
# Creates: find_column()
# ============================================================

find_column <- function(
    column_names,
    exact_candidates,
    regex,
    description
) {
  
  # First prefer an exact column-name match
  exact.match <- exact_candidates[
    exact_candidates %in% column_names
  ]
  
  if (length(exact.match) > 0) {
    return(exact.match[1])
  }
  
  # Otherwise try a case-insensitive regular expression
  regex.match <- grep(
    pattern = regex,
    x = column_names,
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(regex.match) == 0) {
    stop(
      "Could not identify the ",
      description,
      " column.\nAvailable columns are:\n",
      paste(column_names, collapse = ", ")
    )
  }
  
  if (length(regex.match) > 1) {
    message(
      "Multiple possible ",
      description,
      " columns found: ",
      paste(regex.match, collapse = ", "),
      ". Using: ",
      regex.match[1]
    )
  }
  
  regex.match[1]
}

# ============================================================
# 7. Validate the FST column structure now
# ============================================================

for (i in seq_along(comp)) {
  
  current.names <- names(comp[[i]])
  
  fst.column <- find_column(
    column_names = current.names,
    exact_candidates = c(
      "FstWC",
      "fstWC",
      "fst",
      "Fst"
    ),
    regex = "FstWC|(^|[_.])Fst($|[_.])",
    description = "Fst"
  )
  
  sites.column <- find_column(
    column_names = current.names,
    exact_candidates = c(
      "sites",
      "sitesUsed",
      "nSites",
      "numSites",
      "sitesCalled"
    ),
    regex = "sites|nSites|numSites",
    description = "site-count"
  )
  
  if (!"mid" %in% current.names) {
    stop(
      "No 'mid' coordinate column was found for ",
      comparison_labels[i],
      ". Available columns are: ",
      paste(current.names, collapse = ", ")
    )
  }
  
  message(
    comparison_labels[i],
    ": FST = '",
    fst.column,
    "', sites = '",
    sites.column,
    "', position = 'mid'."
  )
}

message(
  "\nInitialization complete. Objects created:\n",
  "  paf.table\n",
  "  best.pairs\n",
  "  comp\n",
  "  comparison_labels\n",
  "  find_column()\n",
  "  chrom_data_chr"
)

# Directory creation
output.dir <- file.path(
  "/Users/mariamadrid/Downloads",
  "FST_synteny_by_chromosome"
)

dir.create(
  output.dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# Convert all four FST comparisons into one long table
# ------------------------------------------------------------

fst.long.list <- vector("list", length(comp))

for (i in seq_along(comp)) {
  
  fst.data <- as.data.frame(comp[[i]])
  
  fst.column <- find_column(
    column_names = names(fst.data),
    exact_candidates = c(
      "FstWC",
      "fstWC",
      "fst",
      "Fst"
    ),
    regex = "FstWC|(^|[_.])Fst($|[_.])",
    description = "Fst"
  )
  
  sites.column <- find_column(
    column_names = names(fst.data),
    exact_candidates = c(
      "sites",
      "sitesUsed",
      "nSites",
      "numSites",
      "sitesCalled"
    ),
    regex = "sites|nSites|numSites",
    description = "site-count"
  )
  
  if (!"mid" %in% names(fst.data)) {
    stop(
      "No 'mid' column was found for ",
      comparison_labels[i]
    )
  }
  
  fst.long.list[[i]] <- fst.data %>%
    transmute(
      chromosome = as.character(scaffold),
      
      # Chromosome-local LW coordinate
      position = suppressWarnings(
        as.numeric(as.character(mid))
      ),
      
      fst = suppressWarnings(
        as.numeric(as.character(.data[[fst.column]]))
      ),
      
      sites = suppressWarnings(
        as.numeric(as.character(.data[[sites.column]]))
      ),
      
      comparison = comparison_labels[i]
    ) %>%
    filter(
      chromosome %in% chr.order,
      !is.na(position),
      !is.na(fst),
      !is.na(sites),
      sites >= 400
    ) %>%
    mutate(
      # Negative FST estimates are conventionally plotted as zero
      fst = pmax(fst, 0)
    )
}

fst.long <- bind_rows(fst.long.list) %>%
  mutate(
    chromosome = factor(
      chromosome,
      levels = chr.order
    ),
    comparison = factor(
      comparison,
      levels = comparison_labels
    )
  )

if (nrow(fst.long) == 0) {
  stop("No FST windows remained after filtering.")
}

# ------------------------------------------------------------
# Check that LW is the PAF target
# ------------------------------------------------------------

if (!all(chr.order %in% unique(paf.table$t.name))) {
  warning(
    "Not every CHR1-CHR11 chromosome was found in paf.table$t.name. ",
    "Check whether the LW assembly is actually the PAF target."
  )
}

# ------------------------------------------------------------
# Optional visual settings
# ------------------------------------------------------------

fst.maximum <- 1
high.fst.threshold <- 0.5

population.colors <- c(
  "Belgium SW/LW"  = "#1B9E77",
  "France SW/LW"   = "#7570B3",
  "Portugal SW/LW" = "#D95F02",
  "Spain SW/LW"    = "#E7298A"
)

# ------------------------------------------------------------
# Generate one combined plot per chromosome
# ------------------------------------------------------------

combined.plots <- vector(
  "list",
  length(chr.order)
)

names(combined.plots) <- chr.order

for (lw.chr in chr.order) {
  
  sw.chr <- best.pairs %>%
    filter(as.character(t.name) == lw.chr) %>%
    pull(q.name) %>%
    as.character()
  
  if (length(sw.chr) != 1) {
    warning(
      "Skipping ", lw.chr,
      ": expected one best SW match, but found ",
      length(sw.chr)
    )
    next
  }
  
  message(
    "Creating ", lw.chr,
    " versus ", sw.chr
  )
  
  # ----------------------------------------------------------
  # Synteny data
  # ----------------------------------------------------------
  
  paf.chr <- paf.table %>%
    filter(
      t.name == lw.chr,
      q.name == sw.chr,
      !is.na(t.start),
      !is.na(t.end)
    )
  
  if (nrow(paf.chr) == 0) {
    warning(
      "Skipping ", lw.chr,
      ": no PAF alignments found for ", sw.chr
    )
    next
  }
  
  # Use the PAF target length as the shared reference range.
  # This is the LW chromosome length when LW is the PAF target.
  chromosome.length <- unique(
    as.numeric(paf.chr$t.len)
  )
  
  chromosome.length <- chromosome.length[
    is.finite(chromosome.length)
  ]
  
  if (length(chromosome.length) == 0) {
    stop(
      "No valid t.len value found for ",
      lw.chr
    )
  }
  
  if (length(unique(chromosome.length)) > 1) {
    warning(
      "Multiple target lengths found for ",
      lw.chr,
      "; using the maximum."
    )
  }
  
  chromosome.length <- max(chromosome.length)
  
  # ----------------------------------------------------------
  # FST panel
  # ----------------------------------------------------------
  
  fst.chr <- fst.long %>%
    filter(
      as.character(chromosome) == lw.chr,
      position >= 0,
      position <= chromosome.length
    )
  
  if (nrow(fst.chr) == 0) {
    warning(
      "No FST data found for ",
      lw.chr
    )
  }
  
  p.fst <- ggplot(
    fst.chr,
    aes(
      x = position,
      y = fst,
      colour = comparison
    )
  ) +
    # Alternating faint window guides
    geom_hline(
      yintercept = c(0.25, 0.5, 0.75),
      linewidth = 0.25,
      colour = "grey85"
    ) +
    
    # Optional threshold for visually high FST
    geom_hline(
      yintercept = high.fst.threshold,
      linewidth = 0.45,
      linetype = "dashed",
      colour = "grey35"
    ) +
    
    geom_point(
      size = 0.55,
      alpha = 0.55
    ) +
    
    facet_grid(
      rows = vars(comparison),
      scales = "fixed",
      switch = "y"
    ) +
    
    scale_colour_manual(
      values = population.colors,
      drop = FALSE
    ) +
    
    scale_x_continuous(
      limits = c(0, chromosome.length),
      expand = expansion(mult = c(0, 0)),
      labels = label_number(
        scale = 1e-6,
        accuracy = 1
      )
    ) +
    
    scale_y_continuous(
      limits = c(0, fst.maximum),
      breaks = c(0, 0.5, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    
    labs(
      x = NULL,
      y = expression(F[ST])
    ) +
    
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      
      # Remove the FST x-axis so there is one coordinate axis
      # at the bottom of the complete figure.
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        angle = 0,
        size = 8,
        hjust = 1
      ),
      
      panel.spacing.y = grid::unit(1.5, "mm"),
      plot.margin = margin(
        t = 5,
        r = 5,
        b = 0,
        l = 5
      )
    )
  
  # ----------------------------------------------------------
  # SVbyEye synteny panel
  # ----------------------------------------------------------
  
  p.synteny <- plotGenome(
    paf.table = paf.chr,
    chromosomes = lw.chr,
    chromosome.bar.width = grid::unit(2, "mm"),
    min.query.aligned.bp = 5000000,
    color.by = "direction",
    color.palette = c(
      "+" = "#3B82F6",
      "-" = "#EF4444"
    )
  ) +
    # Apply exactly the same LW target-coordinate range
    coord_cartesian(
      xlim = c(0, chromosome.length),
      expand = FALSE,
      clip = "on"
    ) +
    labs(
      x = "Position on LW reference chromosome (Mb)"
    ) +
    theme(
      plot.margin = margin(
        t = 0,
        r = 5,
        b = 5,
        l = 5
      )
    )
  
  # ----------------------------------------------------------
  # Stack FST directly over synteny
  # ----------------------------------------------------------
  
  p.combined <- p.fst / p.synteny +
    plot_layout(
      heights = c(2.2, 1.25),
      guides = "collect"
    ) +
    plot_annotation(
      title = paste0(
        "LW ", lw.chr,
        " versus SW ", sw.chr
      ),
      subtitle = paste0(
        "FST and synteny shown in LW ",
        lw.chr,
        " reference coordinates"
      ),
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 14
        ),
        plot.subtitle = element_text(
          size = 10
        )
      )
    )
  
  combined.plots[[lw.chr]] <- p.combined
  
  filename.prefix <- paste0(
    "LW_",
    lw.chr,
    "_vs_SW_",
    sw.chr,
    "_FST_synteny"
  )
  
  ggsave(
    filename = file.path(
      output.dir,
      paste0(filename.prefix, ".pdf")
    ),
    plot = p.combined,
    width = 14,
    height = 9,
    units = "in"
  )
  
  ggsave(
    filename = file.path(
      output.dir,
      paste0(filename.prefix, ".png")
    ),
    plot = p.combined,
    width = 14,
    height = 9,
    units = "in",
    dpi = 300
  )
}

