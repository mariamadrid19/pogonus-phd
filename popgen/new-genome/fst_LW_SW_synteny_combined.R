library(SVbyEye)
library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)

# ============================================================
# Input paths and settings
# ============================================================

paf.file <- "/Users/mariamadrid/Downloads/PchalceusLWrefV2.paf.gz"
lw.fst.directory <- "/Users/mariamadrid/Downloads/LW"
sw.fst.directory <- "/Users/mariamadrid/Downloads/SW-stats"
lw.length.file <- file.path(lw.fst.directory, "scaflengthsLW.txt")
sw.length.file <- file.path(sw.fst.directory, "scaflengthsSW.txt")

chr.order <- paste0("CHR", 1:11)
minimum.sites <- 400
minimum.aligned.bp <- 5000000
fst.maximum <- 1
high.fst.threshold <- 0.5

population.info <- data.frame(
  population = c("Belgium", "France", "Portugal", "Spain"),
  population.number = 0:3,
  comparison = c(
    "Belgium SW/LW",
    "France SW/LW",
    "Portugal SW/LW",
    "Spain SW/LW"
  ),
  stringsAsFactors = FALSE
)

comparison_labels <- population.info$comparison

population.colors <- c(
  "Belgium SW/LW" = "#1B9E77",
  "France SW/LW" = "#7570B3",
  "Portugal SW/LW" = "#D95F02",
  "Spain SW/LW" = "#E7298A"
)

# The first directory contains the LW-Fst + synteny plots.
lw.output.dir <- file.path(
  "/Users/mariamadrid/Downloads",
  "Fst_LW_synteny_by_chromosome"
)

# The second directory contains LW-Fst + synteny + SW-Fst plots.
both.output.dir <- file.path(
  "/Users/mariamadrid/Downloads",
  "Fst_LW_synteny_SW_by_chromosome"
)

dir.create(lw.output.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(both.output.dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 1. Import and validate the PAF alignments
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
  "q.name", "q.len", "q.start", "q.end",
  "strand", "t.name", "t.len", "t.start", "t.end",
  "n.match", "aln.len"
)

missing.paf.columns <- setdiff(required.paf.columns, names(paf.table))

if (length(missing.paf.columns) > 0) {
  stop(
    "Required PAF columns are missing: ",
    paste(missing.paf.columns, collapse = ", ")
  )
}

paf.table <- paf.table %>%
  mutate(
    q.name = as.character(q.name),
    t.name = as.character(t.name)
  )

message(
  "Imported ",
  format(nrow(paf.table), big.mark = ","),
  " PAF alignments."
)

missing.target.chromosomes <- setdiff(
  chr.order,
  unique(paf.table$t.name)
)

if (length(missing.target.chromosomes) > 0) {
  warning(
    "LW chromosomes absent from paf.table$t.name: ",
    paste(missing.target.chromosomes, collapse = ", "),
    ". Confirm that LW is the PAF target/reference."
  )
}

# ============================================================
# 2. Read and validate both assemblies' chromosome lengths
# ============================================================

if (!file.exists(lw.length.file)) {
  stop("LW chromosome-length file does not exist: ", lw.length.file)
}

chrom_data <- read.table(
  lw.length.file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("chrom", "length") %in% names(chrom_data))) {
  stop("The LW length file must contain columns named 'chrom' and 'length'.")
}

chrom_data_chr <- chrom_data %>%
  transmute(
    chrom = as.character(chrom),
    length = suppressWarnings(as.numeric(as.character(length)))
  ) %>%
  filter(chrom %in% chr.order) %>%
  mutate(chromosome = as.integer(sub("^CHR", "", chrom))) %>%
  arrange(chromosome)

if (nrow(chrom_data_chr) != 11) {
  stop(
    "Expected 11 LW chromosomes, but found ",
    nrow(chrom_data_chr),
    ": ",
    paste(chrom_data_chr$chrom, collapse = ", ")
  )
}

if (anyNA(chrom_data_chr$length) || any(chrom_data_chr$length <= 0)) {
  stop("Some LW chromosome lengths are missing or invalid.")
}

if (!file.exists(sw.length.file)) {
  stop("SW chromosome-length file does not exist: ", sw.length.file)
}

sw_chrom_data <- read.table(
  sw.length.file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!all(c("chrom", "length") %in% names(sw_chrom_data))) {
  stop("The SW length file must contain columns named 'chrom' and 'length'.")
}

sw_chrom_data_chr <- sw_chrom_data %>%
  transmute(
    sw.fst.chromosome = as.character(chrom),
    sw.length = suppressWarnings(as.numeric(as.character(length)))
  ) %>%
  filter(sw.fst.chromosome %in% chr.order) %>%
  mutate(
    chromosome.number = as.integer(
      sub("^CHR", "", sw.fst.chromosome)
    )
  ) %>%
  arrange(chromosome.number)

if (nrow(sw_chrom_data_chr) != 11) {
  stop(
    "Expected 11 SW chromosomes, but found ",
    nrow(sw_chrom_data_chr),
    ": ",
    paste(sw_chrom_data_chr$sw.fst.chromosome, collapse = ", ")
  )
}

if (anyNA(sw_chrom_data_chr$sw.length) ||
    any(sw_chrom_data_chr$sw.length <= 0)) {
  stop("Some SW chromosome lengths are missing or invalid.")
}

# Validate the V2 PAF query labels against the SW-stats chromosome labels using
# both chromosome name and exact length. V2 is expected to have corrected names;
# the script stops rather than silently remapping anything if they disagree.
paf.query.lengths <- paf.table %>%
  filter(q.name %in% chr.order) %>%
  distinct(paf.q.name = q.name, paf.q.length = q.len)

if (anyDuplicated(paf.query.lengths$paf.q.name)) {
  stop("A PAF query name is associated with multiple sequence lengths.")
}

paf.to.sw.map <- paf.query.lengths %>%
  left_join(
    sw_chrom_data_chr %>%
      select(sw.fst.chromosome, sw.length),
    by = c("paf.q.name" = "sw.fst.chromosome")
  )

if (anyNA(paf.to.sw.map$sw.length)) {
  stop(
    "These V2 PAF query names are absent from scaflengthsSW.txt: ",
    paste(
      paf.to.sw.map$paf.q.name[is.na(paf.to.sw.map$sw.length)],
      collapse = ", "
    )
  )
}

length.mismatch <- paf.to.sw.map$paf.q.length != paf.to.sw.map$sw.length

if (any(length.mismatch)) {
  stop(
    "V2 PAF query lengths disagree with scaflengthsSW.txt for: ",
    paste(paf.to.sw.map$paf.q.name[length.mismatch], collapse = ", ")
  )
}

paf.to.sw.map <- paf.to.sw.map %>%
  mutate(sw.fst.chromosome = paf.q.name)

message("V2 PAF query names and lengths match the SW FST assembly:")
print(paf.to.sw.map %>% arrange(paf.q.name))

# ============================================================
# 3. Identify the best SW chromosome for every LW chromosome
# ============================================================

pair.summary <- paf.table %>%
  mutate(n.match = suppressWarnings(as.numeric(as.character(n.match)))) %>%
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
  stop("No CHR-to-CHR PAF alignments were found.")
}

best.pairs <- pair.summary %>%
  dplyr::group_by(.data$t.name) %>%
  dplyr::slice_max(
    order_by = .data$matched.bp,
    n = 1,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    chromosome.number = as.integer(
      sub("^CHR", "", .data$t.name)
    ),
    paf.q.name = .data$q.name
  ) %>%
  dplyr::arrange(.data$chromosome.number) %>%
  dplyr::left_join(
    paf.to.sw.map,
    by = "paf.q.name"
  ) %>%
  dplyr::transmute(
    lw.chromosome = .data$t.name,
    paf.q.name = .data$paf.q.name,
    sw.chromosome = .data$sw.fst.chromosome,
    paf.q.length = .data$paf.q.length,
    matched.bp = .data$matched.bp,
    alignment.count = .data$alignment.count
  )

message("Best LW-SW chromosome pairs:")
print(best.pairs)

# ============================================================
# 4. Functions for locating and importing the FST files
# ============================================================

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

locate_stats_files <- function(directory, genome.prefix) {
  if (!dir.exists(directory)) {
    stop("FST directory does not exist: ", directory)
  }
  
  result <- vector("list", nrow(population.info))
  names(result) <- population.info$population
  
  for (i in seq_len(nrow(population.info))) {
    pop.number <- population.info$population.number[i]
    
    pattern <- paste0(
      "^", genome.prefix, "\\.chr_",
      "(1|2|3|4|5|6|7|8|9|10|11)",
      "\\.pop", pop.number, "\\.stats$"
    )
    
    files <- list.files(
      path = directory,
      pattern = pattern,
      full.names = TRUE
    )
    
    files <- sort_stats_files(files)
    
    if (length(files) != 11) {
      stop(
        "Expected 11 files for ", genome.prefix, " ",
        population.info$population[i], ", but found ", length(files),
        ".\nFiles found:\n", paste(files, collapse = "\n")
      )
    }
    
    detected.chromosomes <- as.integer(
      sub(
        ".*\\.chr_([0-9]+)\\.pop[0-9]+\\.stats$",
        "\\1",
        basename(files)
      )
    )
    
    if (!setequal(detected.chromosomes, 1:11)) {
      stop(
        genome.prefix, " ", population.info$population[i],
        " does not contain exactly chromosomes 1-11. Detected: ",
        paste(detected.chromosomes, collapse = ", ")
      )
    }
    
    result[[i]] <- files
  }
  
  result
}

read_stats_files <- function(files) {
  imported <- lapply(files, function(current.file) {
    x <- data.table::fread(current.file, data.table = FALSE)
    x$source_file <- basename(current.file)
    x
  })
  
  result <- data.table::rbindlist(
    imported,
    use.names = TRUE,
    fill = TRUE
  ) %>%
    as.data.frame()
  
  if (nrow(result) == 0) {
    stop("An imported FST dataset contains no rows.")
  }
  
  if (!"scaffold" %in% names(result)) {
    stop(
      "No 'scaffold' column was found. Available columns: ",
      paste(names(result), collapse = ", ")
    )
  }
  
  result$scaffold <- as.character(result$scaffold)
  result
}

# Strictly select ordinary Fst. The pattern deliberately excludes FstWC.
find_fst_column <- function(column.names) {
  matches <- grep(
    pattern = "^Fst($|[_.])",
    x = column.names,
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(matches) == 0) {
    stop(
      "No ordinary Fst column was found. Available columns:\n",
      paste(column.names, collapse = ", ")
    )
  }
  
  if (length(matches) > 1) {
    stop(
      "More than one ordinary Fst column was found: ",
      paste(matches, collapse = ", "),
      ". Select the required column explicitly."
    )
  }
  
  if (grepl("^FstWC", matches[1], ignore.case = TRUE)) {
    stop("Safety check failed: an FstWC column was selected.")
  }
  
  message("Using ordinary Fst column: ", matches[1])
  matches[1]
}

find_sites_column <- function(column.names) {
  exact.candidates <- c(
    "sites", "sitesUsed", "nSites", "numSites", "sitesCalled"
  )
  
  exact.match <- exact.candidates[exact.candidates %in% column.names]
  if (length(exact.match) > 0) return(exact.match[1])
  
  matches <- grep(
    "sites|nSites|numSites",
    column.names,
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(matches) == 0) {
    stop(
      "No site-count column was found. Available columns: ",
      paste(column.names, collapse = ", ")
    )
  }
  
  matches[1]
}

# ============================================================
# 5. Import LW and SW population-statistics files
# ============================================================

lw.file_lists <- locate_stats_files(
  directory = lw.fst.directory,
  genome.prefix = "Pchal_LW"
)

sw.file_lists <- locate_stats_files(
  directory = sw.fst.directory,
  genome.prefix = "Pchal_SW"
)

comp <- lapply(lw.file_lists, read_stats_files)
comp.sw <- lapply(sw.file_lists, read_stats_files)

names(comp) <- comparison_labels
names(comp.sw) <- comparison_labels

# ============================================================
# 5a. Mandatory audit of the exact ordinary-Fst columns
#
# This audit is intentionally performed before plotting. It records the
# selected column and counts numeric ordinary-Fst estimates. FstWC columns are
# explicitly excluded and cannot enter the plotted data.
# ============================================================

audit_fst_inputs <- function(comparison.list, genome.name) {
  audit.rows <- list()
  row.index <- 1L
  
  for (i in seq_along(comparison.list)) {
    x <- as.data.frame(comparison.list[[i]])
    fst.column <- find_fst_column(names(x))
    sites.column <- find_sites_column(names(x))
    
    # Defensive assertion: reject FstWC even if the helper is changed later.
    if (!grepl("^Fst($|[_.])", fst.column, ignore.case = TRUE) ||
        grepl("^FstWC", fst.column, ignore.case = TRUE)) {
      stop(
        "Internal safety check failed: selected column is not ordinary Fst: ",
        fst.column
      )
    }
    
    fst.values <- suppressWarnings(as.numeric(as.character(x[[fst.column]])))
    site.values <- suppressWarnings(as.numeric(as.character(x[[sites.column]])))
    
    for (chr in chr.order) {
      on.chr <- as.character(x$scaffold) == chr
      numeric.fst <- on.chr & is.finite(fst.values)
      plotted.fst <- numeric.fst & is.finite(site.values) &
        site.values >= minimum.sites
      
      audit.rows[[row.index]] <- data.frame(
        genome = genome.name,
        comparison = comparison_labels[i],
        chromosome = chr,
        selected_column = fst.column,
        total_windows = sum(on.chr, na.rm = TRUE),
        numeric_Fst = sum(numeric.fst, na.rm = TRUE),
        missing_Fst = sum(on.chr, na.rm = TRUE) -
          sum(numeric.fst, na.rm = TRUE),
        plotted_Fst = sum(plotted.fst, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      row.index <- row.index + 1L
    }
  }
  
  bind_rows(audit.rows)
}

fst.audit <- bind_rows(
  audit_fst_inputs(comp, "LW"),
  audit_fst_inputs(comp.sw, "SW")
)

audit.file <- file.path(
  both.output.dir,
  "Fst_columns_and_point_counts.tsv"
)

write.table(
  fst.audit,
  file = audit.file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("\nMandatory ordinary-Fst input audit:")
print(fst.audit, row.names = FALSE)
message("Fst audit saved to: ", audit.file)

# ============================================================
# 6. Convert ordinary Fst data to long tables
# ============================================================

make_fst_long <- function(comparison.list, genome.name) {
  result.list <- vector("list", length(comparison.list))
  
  for (i in seq_along(comparison.list)) {
    fst.data <- as.data.frame(comparison.list[[i]])
    fst.column <- find_fst_column(names(fst.data))
    sites.column <- find_sites_column(names(fst.data))
    
    if (!grepl("^Fst($|[_.])", fst.column, ignore.case = TRUE) ||
        grepl("^FstWC", fst.column, ignore.case = TRUE)) {
      stop(
        "Refusing to plot a non-ordinary-Fst column: ",
        fst.column
      )
    }
    
    if (!"mid" %in% names(fst.data)) {
      stop(
        "No 'mid' coordinate column was found for ", genome.name,
        " ", comparison_labels[i], "."
      )
    }
    
    result.list[[i]] <- fst.data %>%
      transmute(
        chromosome = as.character(scaffold),
        position = suppressWarnings(as.numeric(as.character(mid))),
        fst = suppressWarnings(
          as.numeric(as.character(.data[[fst.column]]))
        ),
        sites = suppressWarnings(
          as.numeric(as.character(.data[[sites.column]]))
        ),
        comparison = comparison_labels[i],
        genome = genome.name
      ) %>%
      filter(
        chromosome %in% chr.order,
        !is.na(position),
        !is.na(fst),
        !is.na(sites),
        sites >= minimum.sites
      ) %>%
      mutate(fst = pmax(fst, 0))
  }
  
  bind_rows(result.list) %>%
    mutate(
      chromosome = factor(chromosome, levels = chr.order),
      comparison = factor(comparison, levels = comparison_labels),
      genome = factor(genome, levels = c("LW", "SW"))
    )
}

fst.long.lw <- make_fst_long(comp, "LW")
fst.long.sw <- make_fst_long(comp.sw, "SW")

if (nrow(fst.long.lw) == 0) {
  stop("No LW Fst windows remained after filtering.")
}

if (nrow(fst.long.sw) == 0) {
  stop("No SW Fst windows remained after filtering.")
}

# ============================================================
# 7. Reusable ordinary-Fst panel
# ============================================================

make_fst_panel <- function(
    data,
    plot.limits,
    genome.label,
    show.x.axis = FALSE,
    reversed = FALSE,
    comparison.order = comparison_labels
) {
  # Set facet order independently for each panel. This allows the LW tracks
  # to run Belgium-to-Spain and the SW tracks to mirror them Spain-to-Belgium.
  data <- data %>%
    dplyr::mutate(
      comparison = factor(
        as.character(.data$comparison),
        levels = comparison.order
      )
    )
  
  axis.title <- if (reversed) {
    paste0(genome.label, " position in displayed/reversed orientation (Mb)")
  } else {
    paste0(genome.label, " chromosome position (Mb)")
  }
  
  p <- ggplot(
    data,
    aes(x = plot.position, y = fst, colour = comparison)
  ) +
    geom_hline(
      yintercept = c(0.25, 0.5, 0.75),
      linewidth = 0.25,
      colour = "grey85"
    ) +
    geom_hline(
      yintercept = high.fst.threshold,
      linewidth = 0.45,
      linetype = "dashed",
      colour = "grey35"
    ) +
    geom_point(size = 0.55, alpha = 0.55) +
    facet_grid(
      rows = vars(comparison),
      scales = "fixed",
      switch = "y"
    ) +
    scale_colour_manual(values = population.colors, drop = FALSE) +
    scale_x_continuous(
      limits = plot.limits,
      expand = expansion(mult = c(0, 0)),
      labels = label_number(scale = 1e-6, accuracy = 1)
    ) +
    scale_y_continuous(
      limits = c(0, fst.maximum),
      breaks = c(0, 0.5, 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste0(genome.label, " Fst"),
      x = if (show.x.axis) axis.title else NULL,
      y = expression(F[ST])
    ) +
    theme_classic(base_size = 10) +
    theme(
      legend.position = "none",
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(angle = 0, size = 8, hjust = 1),
      panel.spacing.y = grid::unit(1.5, "mm"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0),
      plot.margin = margin(t = 2, r = 5, b = 2, l = 5)
    )
  
  if (!show.x.axis) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    )
  }
  
  p
}

# ============================================================
# 8. Exact full-chromosome coordinates and synteny panel
# ============================================================

get_pair_coordinates <- function(paf.chr, lw.length, sw.length) {
  paf.lw.length <- unique(as.numeric(paf.chr$t.len))
  paf.sw.length <- unique(as.numeric(paf.chr$q.len))
  
  if (length(paf.lw.length) != 1 || paf.lw.length != lw.length) {
    stop(
      "PAF target length does not match scaflengthsLW.txt: PAF=",
      paste(paf.lw.length, collapse = ","),
      ", length file=", lw.length
    )
  }
  
  if (length(paf.sw.length) != 1 || paf.sw.length != sw.length) {
    stop(
      "PAF query length does not match its length-mapped SW chromosome: PAF=",
      paste(paf.sw.length, collapse = ","),
      ", scaflengthsSW.txt=", sw.length
    )
  }
  
  forward.bp <- sum(
    paf.chr$aln.len[paf.chr$strand == "+"],
    na.rm = TRUE
  )
  reverse.bp <- sum(
    paf.chr$aln.len[paf.chr$strand == "-"],
    na.rm = TRUE
  )
  
  flip.query <- reverse.bp > forward.bp
  
  list(
    lw.length = lw.length,
    sw.length = sw.length,
    flip.query = flip.query,
    plot.limits = c(0, max(lw.length, sw.length, na.rm = TRUE))
  )
}

# Draw synteny directly from raw PAF coordinates. Unlike plotGenome(), this
# draws chromosome bars from the authoritative sequence lengths rather than
# trimming them to the min/max aligned positions. Therefore FST and synteny
# use exactly the same 0-to-chromosome-length coordinate system.
make_exact_synteny_panel <- function(
    paf.chr,
    coords,
    lw.chr,
    sw.chr,
    show.x.axis = FALSE
) {
  ribbons <- paf.chr %>%
    dplyr::mutate(
      alignment.id = dplyr::row_number(),
      
      # Apply the same whole-query orientation choice used previously.
      q.start.display = if (coords$flip.query) {
        coords$sw.length - .data$q.end
      } else {
        .data$q.start
      },
      q.end.display = if (coords$flip.query) {
        coords$sw.length - .data$q.start
      } else {
        .data$q.end
      },
      display.strand = if (coords$flip.query) {
        ifelse(.data$strand == "+", "-", "+")
      } else {
        as.character(.data$strand)
      },
      
      # Query coordinates corresponding to the target start and target end.
      q.at.target.start = ifelse(
        .data$display.strand == "+",
        .data$q.start.display,
        .data$q.end.display
      ),
      q.at.target.end = ifelse(
        .data$display.strand == "+",
        .data$q.end.display,
        .data$q.start.display
      )
    )
  
  # Cubic Bezier edge with vertical tangents at both chromosome bars. This
  # recreates the smooth S-shaped ribbon boundaries used by SVbyEye while
  # retaining the unshifted, full-chromosome coordinate system.
  bezier_edge <- function(x.start, y.start, x.end, y.end, n = 40) {
    t <- seq(0, 1, length.out = n)
    
    # Control points retain the x coordinate of their nearest endpoint.
    # This creates a smooth departure from and arrival at each chromosome bar.
    control1.x <- x.start
    control1.y <- y.start + (y.end - y.start) * 0.35
    control2.x <- x.end
    control2.y <- y.start + (y.end - y.start) * 0.65
    
    data.frame(
      x =
        (1 - t)^3 * x.start +
        3 * (1 - t)^2 * t * control1.x +
        3 * (1 - t) * t^2 * control2.x +
        t^3 * x.end,
      y =
        (1 - t)^3 * y.start +
        3 * (1 - t)^2 * t * control1.y +
        3 * (1 - t) * t^2 * control2.y +
        t^3 * y.end
    )
  }
  
  ribbon.polygons <- dplyr::bind_rows(lapply(
    seq_len(nrow(ribbons)),
    function(i) {
      # Right boundary: target end down to its corresponding query endpoint.
      right.edge <- bezier_edge(
        x.start = ribbons$t.end[i],
        y.start = 2,
        x.end = ribbons$q.at.target.end[i],
        y.end = 1
      )
      
      # Left boundary: corresponding query endpoint back up to target start.
      left.edge <- bezier_edge(
        x.start = ribbons$q.at.target.start[i],
        y.start = 1,
        x.end = ribbons$t.start[i],
        y.end = 2
      )
      
      polygon <- dplyr::bind_rows(
        data.frame(x = ribbons$t.start[i], y = 2),
        right.edge,
        data.frame(x = ribbons$q.at.target.start[i], y = 1),
        left.edge
      )
      
      polygon$alignment.id <- ribbons$alignment.id[i]
      polygon$direction <- ribbons$display.strand[i]
      polygon$vertex.order <- seq_len(nrow(polygon))
      polygon
    }
  ))
  
  chromosome.bars <- data.frame(
    genome = c(paste0("LW ", lw.chr), paste0("SW ", sw.chr)),
    xmin = c(0, 0),
    xmax = c(coords$lw.length, coords$sw.length),
    y = c(2, 1),
    fill = c("#B45309", "grey75")
  )
  
  p <- ggplot() +
    geom_polygon(
      data = ribbon.polygons,
      aes(
        x = .data$x,
        y = .data$y,
        group = .data$alignment.id,
        order = .data$vertex.order,
        fill = .data$direction
      ),
      alpha = 0.48,
      colour = NA
    ) +
    scale_fill_manual(
      values = c("+" = "#3B82F6", "-" = "#EF4444"),
      name = "Alignment\ndirection"
    ) +
    ggnewscale::new_scale_fill() +
    # Use SVbyEye's rounded chromosome geometry, but with a much thinner
    # fixed physical width. The thickness therefore stays consistent when
    # the complete multi-panel figure is resized.
    SVbyEye:::geom_roundrect(
      data = chromosome.bars,
      aes(
        xmin = .data$xmin,
        xmax = .data$xmax,
        y = .data$y,
        fill = .data$genome
      ),
      rect_height = grid::unit(1.2, "mm"),
      radius = grid::unit(0.6, "mm"),
      colour = "black",
      linewidth = 0.3,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = setNames(chromosome.bars$fill, chromosome.bars$genome),
      name = "Chromosome"
    ) +
    scale_x_continuous(
      limits = coords$plot.limits,
      expand = expansion(mult = c(0, 0)),
      labels = label_number(scale = 1e-6, accuracy = 1)
    ) +
    scale_y_continuous(
      limits = c(0.94, 2.06),
      breaks = c(2, 1),
      labels = c(paste0("LW ", lw.chr), paste0("SW ", sw.chr)),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = if (show.x.axis) "Chromosome position (Mb)" else NULL,
      y = NULL
    ) +
    theme_classic(base_size = 10) +
    theme(
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 9),
      legend.position = "right",
      plot.margin = margin(t = 1, r = 5, b = 1, l = 5)
    )
  
  if (!show.x.axis) {
    p <- p + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank()
    )
  }
  
  p
}

# ============================================================
# PART A: Original layout -- LW Fst above synteny
# ============================================================

lw.synteny.plots <- setNames(vector("list", length(chr.order)), chr.order)

for (lw.chr in chr.order) {
  pair.row <- best.pairs %>% filter(lw.chromosome == lw.chr)
  
  paf.q.chr <- as.character(pair.row$paf.q.name)
  sw.chr <- as.character(pair.row$sw.chromosome)
  
  if (nrow(pair.row) != 1 || length(paf.q.chr) != 1 || length(sw.chr) != 1) {
    warning("Skipping ", lw.chr, ": expected one SW match.")
    next
  }
  
  paf.chr <- paf.table %>%
    filter(
      t.name == lw.chr,
      q.name == paf.q.chr,
      !is.na(t.start), !is.na(t.end),
      !is.na(q.start), !is.na(q.end)
    )
  
  if (nrow(paf.chr) == 0) {
    warning("Skipping ", lw.chr, ": no PAF alignments for ", paf.q.chr)
    next
  }
  
  lw.authoritative.length <- chrom_data_chr$length[
    chrom_data_chr$chrom == lw.chr
  ]
  sw.authoritative.length <- sw_chrom_data_chr$sw.length[
    sw_chrom_data_chr$sw.fst.chromosome == sw.chr
  ]
  
  coords <- get_pair_coordinates(
    paf.chr,
    lw.length = lw.authoritative.length,
    sw.length = sw.authoritative.length
  )
  
  fst.chr.lw <- fst.long.lw %>%
    filter(
      as.character(chromosome) == lw.chr,
      position >= 0,
      position <= coords$lw.length
    ) %>%
    mutate(
      plot.position = position
    ) %>%
    filter(
      plot.position >= coords$plot.limits[1],
      plot.position <= coords$plot.limits[2]
    )
  
  p.fst.lw <- make_fst_panel(
    data = fst.chr.lw,
    plot.limits = coords$plot.limits,
    genome.label = paste0("LW ", lw.chr),
    show.x.axis = FALSE
  )
  
  p.synteny <- make_exact_synteny_panel(
    paf.chr = paf.chr,
    coords = coords,
    lw.chr = lw.chr,
    sw.chr = sw.chr,
    show.x.axis = TRUE
  )
  
  p.combined <- p.fst.lw / p.synteny +
    plot_layout(heights = c(2.2, 1.25), guides = "collect") +
    plot_annotation(
      title = paste0("LW ", lw.chr, " versus SW ", sw.chr),
      subtitle = "LW FST above chromosome synteny",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10)
      )
    )
  
  lw.synteny.plots[[lw.chr]] <- p.combined
  
  prefix <- paste0(
    "LW_", lw.chr, "_vs_SW_", sw.chr, "_Fst_synteny"
  )
  
  ggsave(
    file.path(lw.output.dir, paste0(prefix, ".pdf")),
    p.combined,
    width = 14,
    height = 9,
    units = "in"
  )
  
  ggsave(
    file.path(lw.output.dir, paste0(prefix, ".png")),
    p.combined,
    width = 14,
    height = 9,
    units = "in",
    dpi = 300
  )
}

# Multipage PDF for Part A
pdf(
  file.path(lw.output.dir, "all_chromosomes_LW_Fst_synteny.pdf"),
  width = 14,
  height = 9,
  onefile = TRUE
)
for (lw.chr in chr.order) {
  if (!is.null(lw.synteny.plots[[lw.chr]])) {
    print(lw.synteny.plots[[lw.chr]])
  }
}
dev.off()

# ============================================================
# PART B: LW Fst above, synteny middle, SW Fst below
# ============================================================

both.fst.plots <- setNames(vector("list", length(chr.order)), chr.order)

for (lw.chr in chr.order) {
  pair.row <- best.pairs %>% filter(lw.chromosome == lw.chr)
  
  paf.q.chr <- as.character(pair.row$paf.q.name)
  sw.chr <- as.character(pair.row$sw.chromosome)
  
  if (nrow(pair.row) != 1 || length(paf.q.chr) != 1 || length(sw.chr) != 1) {
    warning("Skipping ", lw.chr, ": expected one SW match.")
    next
  }
  
  message(
    "Creating LW ", lw.chr,
    " versus SW ", sw.chr,
    " (PAF query label ", paf.q.chr, ")"
  )
  
  paf.chr <- paf.table %>%
    filter(
      t.name == lw.chr,
      q.name == paf.q.chr,
      !is.na(t.start), !is.na(t.end),
      !is.na(q.start), !is.na(q.end)
    )
  
  if (nrow(paf.chr) == 0) {
    warning("Skipping ", lw.chr, ": no PAF alignments for ", paf.q.chr)
    next
  }
  
  lw.authoritative.length <- chrom_data_chr$length[
    chrom_data_chr$chrom == lw.chr
  ]
  sw.authoritative.length <- sw_chrom_data_chr$sw.length[
    sw_chrom_data_chr$sw.fst.chromosome == sw.chr
  ]
  
  coords <- get_pair_coordinates(
    paf.chr,
    lw.length = lw.authoritative.length,
    sw.length = sw.authoritative.length
  )
  
  message(
    "  SW display orientation: ",
    ifelse(coords$flip.query, "reversed", "forward")
  )
  
  fst.chr.lw <- fst.long.lw %>%
    filter(
      as.character(chromosome) == lw.chr,
      position >= 0,
      position <= coords$lw.length
    ) %>%
    mutate(
      plot.position = position
    ) %>%
    filter(
      plot.position >= coords$plot.limits[1],
      plot.position <= coords$plot.limits[2]
    )
  
  fst.chr.sw <- fst.long.sw %>%
    filter(
      as.character(chromosome) == sw.chr,
      position >= 0,
      position <= coords$sw.length
    ) %>%
    mutate(
      plot.position = if (coords$flip.query) {
        coords$sw.length - position
      } else {
        position
      }
    ) %>%
    filter(
      plot.position >= coords$plot.limits[1],
      plot.position <= coords$plot.limits[2]
    )
  
  if (nrow(fst.chr.lw) == 0) {
    warning("No LW Fst data found for ", lw.chr)
  }
  if (nrow(fst.chr.sw) == 0) {
    warning("No SW Fst data found for ", sw.chr)
  }
  
  p.fst.lw <- make_fst_panel(
    data = fst.chr.lw,
    plot.limits = coords$plot.limits,
    genome.label = paste0("LW ", lw.chr),
    show.x.axis = FALSE
  )
  
  p.synteny <- make_exact_synteny_panel(
    paf.chr = paf.chr,
    coords = coords,
    lw.chr = lw.chr,
    sw.chr = sw.chr,
    show.x.axis = FALSE
  )
  
  p.fst.sw <- make_fst_panel(
    data = fst.chr.sw,
    plot.limits = coords$plot.limits,
    genome.label = paste0("SW ", sw.chr),
    show.x.axis = TRUE,
    reversed = coords$flip.query,
    comparison.order = rev(comparison_labels)
  )
  
  p.combined <- p.fst.lw / p.synteny / p.fst.sw +
    plot_layout(
      heights = c(2.2, 1.25, 2.2),
      guides = "collect"
    ) +
    plot_annotation(
      title = paste0("LW ", lw.chr, " versus SW ", sw.chr),
      subtitle = if (coords$flip.query) {
        paste0(
          "FST for LW above and SW below; SW ",
          sw.chr, " is displayed in reverse orientation"
        )
      } else {
        "FST for LW above and SW below"
      },
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10)
      )
    )
  
  both.fst.plots[[lw.chr]] <- p.combined
  
  prefix <- paste0(
    "LW_", lw.chr, "_vs_SW_", sw.chr,
    "_LW_Fst_synteny_SW_Fst"
  )
  
  ggsave(
    file.path(both.output.dir, paste0(prefix, ".pdf")),
    p.combined,
    width = 14,
    height = 14,
    units = "in"
  )
  
  ggsave(
    file.path(both.output.dir, paste0(prefix, ".png")),
    p.combined,
    width = 14,
    height = 14,
    units = "in",
    dpi = 300
  )
}

# Multipage PDF for Part B
pdf(
  file.path(
    both.output.dir,
    "all_chromosomes_LW_Fst_synteny_SW_Fst.pdf"
  ),
  width = 14,
  height = 14,
  onefile = TRUE
)
for (lw.chr in chr.order) {
  if (!is.null(both.fst.plots[[lw.chr]])) {
    print(both.fst.plots[[lw.chr]])
  }
}
dev.off()

message("Finished.")
message("Part A output: ", lw.output.dir)
message("Part B output: ", both.output.dir)
