library(SVbyEye)
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(scales)
library(ggnewscale)

# ============================================================
# Inputs
# PAF target = LW; PAF query = SW
# ============================================================

root.dir <- "/Users/mariamadrid/Documents/phd-pogonus"
synteny.dir <- file.path(root.dir, "synteny-fst-repeats")
fst.dirs <- c(
  LW = file.path(root.dir, "popgen/LW-stats-2026"),
  SW = file.path(root.dir, "popgen/SW-stats-2026")
)
repeat.dir <- file.path(root.dir, "repeat-analysis-2026")
expression.dirs <- c(
  LW = "/Users/mariamadrid/Downloads/enrichment-LW-2026",
  SW = "/Users/mariamadrid/Downloads/enrichment-SW-2026"
)
genetic.map.dir <- "/Users/mariamadrid/Downloads/all-recombination-maps"

paf.file <- file.path(synteny.dir, "PchalceusLWrefV2.paf.gz")
lw.agp.file <- file.path(synteny.dir, "Pchalceus_LW_final.fasta.agp")
sw.agp.file <- file.path(synteny.dir, "Pchalceus_SW_final.fasta.agp")
repeat.files <- c(
  LW = file.path(repeat.dir, "PchalceusLW_normalized_repeat_content.tsv"),
  SW = file.path(repeat.dir, "PchalceusSW_normalized_repeat_content.tsv")
)
gtf.files <- c(
  LW = file.path(expression.dirs[["LW"]], "braker_LW.gtf"),
  SW = file.path(expression.dirs[["SW"]], "braker_SW.gtf")
)
de.files <- c(
  LW = file.path(expression.dirs[["LW"]], "DESeq2_SW_vs_LW.tsv"),
  SW = file.path(expression.dirs[["SW"]], "DESeq2_SW_vs_LW.tsv")
)
genetic.map.file <- file.path(
  genetic.map.dir,
  "all_recombination_rates_cM_per_Mb.tsv"
)

mask.high.repeats <- TRUE
repeat.mask.threshold <- 0.75
repeat.mask.fill <- "white"
repeat.mask.alpha <- 0.80

output.dir <- file.path(
  synteny.dir,
  "LW_SW_FST_repeats_genes_AGP_synteny_75pct_white_mask"
)
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

chr.order <- paste0("CHR", 1:11)
flip.sw.agp.chromosomes <- c("CHR3", "CHR6", "CHR11")
fst.maximum <- 1
high.fst.threshold <- 0.5
de.fdr <- 0.05
de.lfc <- 0

population.info <- data.frame(
  population = c("Belgium", "France", "Portugal", "Spain"),
  pop = 0:3,
  comparison = c("Belgium SW/LW", "France SW/LW",
                 "Portugal SW/LW", "Spain SW/LW")
)
comparison.labels <- population.info$comparison
population.colors <- c(
  "Belgium SW/LW" = "#1B9E77",
  "France SW/LW" = "#7570B3",
  "Portugal SW/LW" = "#D95F02",
  "Spain SW/LW" = "#E7298A"
)
direction.colors <- c(
  "SW up" = "#0072B2",
  "LW up" = "#FF0054",
  "Not significant" = "grey70"
)

required.files <- c(paf.file, lw.agp.file, sw.agp.file,
                    repeat.files, gtf.files, de.files, genetic.map.file)
missing.files <- required.files[!file.exists(required.files)]
if (length(missing.files)) {
  stop("Missing input files:\n", paste(missing.files, collapse = "\n"))
}

# ============================================================
# Helpers
# ============================================================

read.agp <- function(file) {
  x <- read.table(file, header = FALSE, sep = "\t", quote = "",
                  comment.char = "#", fill = TRUE,
                  stringsAsFactors = FALSE)
  if (ncol(x) < 9) stop("AGP has fewer than 9 columns: ", file)
  x <- x[, 1:9]
  names(x) <- c("chrom", "start", "end", "part", "type",
                "component", "component.start", "component.end",
                "orientation")
  x %>%
    mutate(chrom = toupper(as.character(chrom)),
           start = as.numeric(start), end = as.numeric(end)) %>%
    filter(type == "W", chrom %in% chr.order,
           is.finite(start), is.finite(end)) %>%
    arrange(match(chrom, chr.order), start) %>%
    group_by(chrom) %>%
    mutate(scaffold.number = row_number(),
           scaffold.color = if_else(scaffold.number %% 2 == 1,
                                    "dark_blue", "light_blue")) %>%
    ungroup()
}

read.fst <- function(directory, assembly) {
  out <- vector("list", nrow(population.info))
  for (i in seq_len(nrow(population.info))) {
    pop <- population.info$pop[i]
    pattern <- paste0("^Pchal_", assembly,
                      "\\.chr_(1|2|3|4|5|6|7|8|9|10|11)",
                      "\\.pop", pop, "\\.stats$")
    files <- list.files(directory, pattern = pattern, full.names = TRUE)
    number <- as.integer(sub(".*\\.chr_([0-9]+)\\.pop[0-9]+\\.stats$",
                             "\\1", basename(files)))
    files <- files[order(number)]
    if (length(files) != 11 || !setequal(number, 1:11)) {
      stop("Expected 11 FST files for pop", pop, "; found ", length(files))
    }
    x <- as.data.frame(rbindlist(lapply(files, fread),
                                 use.names = TRUE, fill = TRUE))
    fst.col <- grep("^Fst(?:$|_)", names(x), value = TRUE,
                    ignore.case = TRUE, perl = TRUE)
    fst.col <- fst.col[!grepl("^FstWC", fst.col, ignore.case = TRUE)][1]
    sites.col <- grep("sites|nSites|numSites", names(x), value = TRUE,
                      ignore.case = TRUE)[1]
    if (is.na(fst.col) || is.na(sites.col) ||
        !all(c("scaffold", "mid") %in% names(x))) {
      stop("Could not identify FST, sites, scaffold, or mid columns for pop", pop)
    }
    message(assembly, " ", population.info$comparison[i], ": using ", fst.col)
    out[[i]] <- x %>%
      transmute(chrom = as.character(scaffold),
                position = as.numeric(mid), fst = as.numeric(.data[[fst.col]]),
                sites = as.numeric(.data[[sites.col]]),
                comparison = population.info$comparison[i]) %>%
      filter(chrom %in% chr.order, is.finite(position), is.finite(fst),
             is.finite(sites), sites >= 400) %>%
      mutate(fst = pmax(fst, 0))
  }
  bind_rows(out) %>%
    mutate(comparison = factor(comparison, levels = comparison.labels))
}

make.fst.plot <- function(x, chromosome.length, assembly,
                          show.x = FALSE,
                          comparison.order = comparison.labels) {
  x <- x %>%
    mutate(comparison = factor(as.character(comparison),
                               levels = comparison.order))
  ggplot(x, aes(position, fst, colour = comparison)) +
    geom_hline(yintercept = c(0.25, 0.5, 0.75), linewidth = 0.25,
               colour = "grey85") +
    geom_hline(yintercept = high.fst.threshold, linewidth = 0.45,
               linetype = "dashed", colour = "grey35") +
    geom_point(aes(alpha = fst), size = 0.55) +
    facet_grid(rows = vars(comparison), switch = "y") +
    scale_colour_manual(values = population.colors, drop = FALSE) +
    scale_alpha_continuous(
      range = c(0.15, 0.9),
      limits = c(0, fst.maximum),
      guide = "none"
    ) +
    scale_x_continuous(limits = c(0, chromosome.length),
                       expand = expansion(mult = c(0, 0)),
                       labels = label_number(scale = 1e-6, accuracy = 1)) +
    scale_y_continuous(limits = c(0, fst.maximum), breaks = c(0, 0.5, 1),
                       expand = expansion(mult = c(0, 0.02))) +
    labs(x = if (show.x) paste0("Position on ", assembly,
                                " reference chromosome (Mb)") else NULL,
         y = expression(F[ST])) +
    theme_classic(base_size = 10) +
    theme(legend.position = "none",
          axis.text.x = if (show.x) element_text() else element_blank(),
          axis.ticks.x = if (show.x) element_line() else element_blank(),
          axis.line.x = if (show.x) element_line() else element_blank(),
          strip.placement = "outside", strip.background = element_blank(),
          strip.text.y.left = element_text(angle = 0, size = 8, hjust = 1),
          panel.spacing.y = grid::unit(1.5, "mm"),
          plot.margin = margin(2, 5, 0, 5))
}

make.repeat.plot <- function(x, chromosome.length, assembly) {
  ggplot(x, aes(position, normalized_coverage)) +
    geom_point(colour = if (assembly == "LW") "#CD9600" else "#6A3D9A",
               size = 0.55, alpha = 0.55) +
    scale_x_continuous(limits = c(0, chromosome.length),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 1),
                       labels = percent_format(accuracy = 1),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = paste0(assembly, " repeats")) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), axis.text.y = element_text(size = 7),
          axis.title.y = element_text(size = 8),
          plot.margin = margin(1, 5, 1, 5))
}

make.genetic.map.plot <- function(x, chromosome.length, assembly) {
  hybrid.color <- if (assembly == "LW") "#FF0054" else "#0072B2"
  
  ggplot() +
    geom_line(
      data = filter(x, family_type == "pure"),
      aes(x = map.position, y = map.value, group = dataset_id),
      colour = "grey75",
      linewidth = 0.4,
      alpha = 0.75
    ) +
    geom_point(
      data = filter(x, family_type == "pure"),
      aes(x = map.position, y = map.value),
      colour = "grey75",
      size = 0.5,
      alpha = 0.65
    ) +
    geom_line(
      data = filter(x, family_type == "hybrid"),
      aes(x = map.position, y = map.value, group = dataset_id),
      colour = hybrid.color,
      linewidth = 0.55,
      alpha = 0.85
    ) +
    geom_point(
      data = filter(x, family_type == "hybrid"),
      aes(x = map.position, y = map.value),
      colour = hybrid.color,
      size = 0.55,
      alpha = 0.8
    ) +
    scale_x_continuous(
      limits = c(0, chromosome.length),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.04))
    ) +
    labs(x = NULL, y = paste0(assembly, " recomb.\n(cM/Mb)")) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      plot.margin = margin(1, 5, 1, 5)
    )
}

merge.mask.intervals <- function(x) {
  if (!nrow(x)) return(data.frame(xmin = numeric(), xmax = numeric()))
  x <- x %>%
    transmute(xmin = pmax(as.numeric(start) - 1, 0),
              xmax = as.numeric(end)) %>%
    filter(is.finite(xmin), is.finite(xmax), xmax >= xmin) %>%
    arrange(xmin, xmax)
  if (!nrow(x)) return(data.frame(xmin = numeric(), xmax = numeric()))
  previous.maximum <- dplyr::lag(cummax(x$xmax), default = -Inf)
  x %>%
    mutate(mask.group = cumsum(xmin > previous.maximum)) %>%
    group_by(mask.group) %>%
    summarise(xmin = min(xmin), xmax = max(xmax), .groups = "drop") %>%
    dplyr::select(xmin, xmax)
}

add.repeat.mask <- function(plot, mask.intervals,
                            ymin = -Inf, ymax = Inf) {
  if (!mask.high.repeats || !nrow(mask.intervals)) return(plot)
  plot +
    geom_rect(
      data = mask.intervals,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = repeat.mask.fill,
      alpha = repeat.mask.alpha,
      colour = NA
    )
}

make.gene.plot <- function(x, chromosome.length, assembly) {
  # Show every gene as a short vertical grey line.
  ggplot() +
    geom_segment(
      data = x,
      aes(x = (start + end) / 2, xend = (start + end) / 2,
          y = 0.375, yend = 0.625),
      colour = "grey65",
      linewidth = 0.32,
      alpha = 0.8
    ) +
    scale_x_continuous(limits = c(0, chromosome.length),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0, 1), breaks = NULL,
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = paste0(assembly, " genes")) +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.line.x = element_blank(), axis.title.y = element_text(size = 8),
          plot.margin = margin(1, 5, 1, 5))
}

make.deg.plot <- function(x, chromosome.length) {
  # Only significant DEGs are shown. Point area represents significance:
  # a smaller raw p-value gives a larger -log10(p-value) and a larger point.
  deg <- x %>%
    filter(direction != "Not significant",
           is.finite(log2FoldChange), is.finite(pvalue), pvalue >= 0) %>%
    mutate(neglog10.pvalue = -log10(pmax(pvalue, .Machine$double.xmin)))
  
  ggplot(deg, aes(x = (start + end) / 2, y = log2FoldChange)) +
    geom_hline(yintercept = 0, linewidth = 0.35,
               linetype = "dashed", colour = "grey45") +
    geom_point(
      aes(colour = direction, size = neglog10.pvalue),
      alpha = 0.72
    ) +
    scale_colour_manual(
      values = direction.colors[c("SW up", "LW up")],
      name = "Differential\nexpression",
      drop = FALSE
    ) +
    scale_size_continuous(
      range = c(0.7, 4.2),
      name = expression(-log[10](italic(p)))
    ) +
    scale_x_continuous(limits = c(0, chromosome.length),
                       expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = expression(log[2](FC))) +
    theme_classic(base_size = 9) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(size = 7),
      axis.title.y = element_text(size = 8),
      plot.margin = margin(1, 5, 1, 5)
    )
}

make.synteny.plot <- function(paf.chr, lw.blocks, sw.blocks,
                              lw.chr, sw.chr, lw.length, sw.length) {
  # Use SVbyEye itself so ribbon geometry and PAF-row retention are exactly
  # the same as in the reference plotting script. min.query.aligned.bp is a
  # query-level threshold; it is NOT a per-alignment length filter.
  p <- SVbyEye::plotGenome(
    paf.table = paf.chr,
    chromosomes = lw.chr,
    chromosome.bar.width = grid::unit(2, "mm"),
    min.query.aligned.bp = 5000000,
    color.by = "direction",
    color.palette = c(
      "+" = "#9CA3AF",
      "-" = "#EF4444"
    )
  )
  
  # Retain SVbyEye's ribbons exactly as drawn, but remove its solid
  # chromosome-bar layers. The AGP components below become the only bars.
  is.svbyeye.chromosome.bar <- vapply(
    p$layers,
    function(layer) inherits(layer$geom, "GeomRoundRect"),
    logical(1)
  )
  p$layers <- p$layers[!is.svbyeye.chromosome.bar]
  
  # These SW chromosomes are reversed relative to the AGP coordinate order.
  # Flip only their SW AGP component coordinates; LW and PAF stay unchanged.
  if (sw.chr %in% flip.sw.agp.chromosomes) {
    sw.blocks <- sw.blocks %>%
      mutate(
        old.start = start,
        old.end = end,
        start = sw.length - old.end + 1,
        end = sw.length - old.start + 1,
        orientation = case_when(
          orientation == "+" ~ "-",
          orientation == "-" ~ "+",
          TRUE ~ orientation
        )
      ) %>%
      arrange(start) %>%
      mutate(
        scaffold.number = row_number(),
        scaffold.color = if_else(
          scaffold.number %% 2 == 1,
          "dark_blue",
          "light_blue"
        )
      ) %>%
      dplyr::select(-old.start, -old.end)
  }
  
  # Overlay alternating AGP components on the rounded chromosome bars.
  # plotGenome places the LW target at y=2 and the SW query at y=1.
  blocks <- bind_rows(
    lw.blocks %>%
      transmute(component, xmin = start - 1, xmax = end,
                y = 2, scaffold.color),
    sw.blocks %>%
      transmute(component, xmin = start - 1, xmax = end,
                y = 1, scaffold.color)
  )
  
  p +
    ggnewscale::new_scale_fill() +
    SVbyEye:::geom_roundrect(
      data = blocks,
      aes(xmin = xmin, xmax = xmax, y = y, fill = scaffold.color),
      rect_height = grid::unit(2.8, "mm"),
      radius = grid::unit(1.4, "mm"),
      colour = "white",
      linewidth = 0.22,
      inherit.aes = FALSE
    ) +
    scale_fill_manual(
      values = c(dark_blue = "#525252", light_blue = "#BDBDBD"),
      guide = "none"
    ) +
    coord_cartesian(
      xlim = c(0, max(lw.length, sw.length)),
      ylim = c(0.88, 2.12),
      expand = FALSE,
      clip = "off"
    ) +
    labs(x = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      plot.margin = margin(0, 5, 0, 5)
    )
}

# ============================================================
# Import data
# ============================================================

paf <- SVbyEye::readPaf(paf.file, include.paf.tags = FALSE)
required.paf <- c("q.name", "q.len", "q.start", "q.end", "strand",
                  "t.name", "t.len", "t.start", "t.end", "n.match", "aln.len")
if (length(setdiff(required.paf, names(paf)))) stop("Required PAF columns missing")

best.pairs <- paf %>%
  filter(t.name %in% chr.order, q.name %in% chr.order) %>%
  group_by(t.name, q.name) %>%
  summarise(matched.bp = sum(as.numeric(n.match), na.rm = TRUE), .groups = "drop") %>%
  group_by(t.name) %>% slice_max(matched.bp, n = 1, with_ties = FALSE) %>%
  ungroup()

fst <- list(
  LW = read.fst(fst.dirs[["LW"]], "LW"),
  SW = read.fst(fst.dirs[["SW"]], "SW")
)

read.repeats <- function(file) {
  fread(file, data.table = FALSE) %>%
    mutate(seqnames = as.character(seqnames), start = as.numeric(start),
           end = as.numeric(end), normalized_coverage = as.numeric(normalized_coverage),
           position = if ("mid" %in% names(.)) as.numeric(mid) else (start + end) / 2) %>%
    filter(seqnames %in% chr.order, is.finite(position),
           is.finite(normalized_coverage))
}
repeats <- list(
  LW = read.repeats(repeat.files[["LW"]]),
  SW = read.repeats(repeat.files[["SW"]])
)
genetic.maps <- fread(genetic.map.file, data.table = FALSE) %>%
  transmute(
    rate_version = as.character(rate_version),
    dataset_id = as.character(dataset_id),
    family_type = as.character(family_type),
    reference_genome = as.character(reference_genome),
    chrom = as.character(chrom),
    map.position = as.numeric(window_midpoint_bp),
    map.value = as.numeric(rate_cM_Mb),
    sufficient_data = as.logical(sufficient_data)
  ) %>%
  filter(
    rate_version == "pure_center_points_hybrid_all_points",
    reference_genome %in% c("LW", "SW"),
    family_type %in% c("pure", "hybrid"),
    chrom %in% chr.order,
    is.finite(map.position),
    is.finite(map.value),
    sufficient_data
  ) %>%
  arrange(reference_genome, chrom, dataset_id, map.position)
agp.lw <- read.agp(lw.agp.file)
agp.sw <- read.agp(sw.agp.file)

read.genes.and.de <- function(gtf.file, de.file, assembly) {
  gtf <- fread(gtf.file, sep = "\t", header = FALSE,
               data.table = FALSE, select = c(1, 3, 4, 5, 7, 9))
  names(gtf) <- c("chrom", "feature", "start", "end", "strand", "attribute")
  genes <- gtf %>%
    filter(feature == "gene", chrom %in% chr.order) %>%
    transmute(chrom, start = as.numeric(start), end = as.numeric(end), strand,
              gene_id = str_remove(str_trim(attribute), ";$")) %>%
    distinct(gene_id, .keep_all = TRUE)
  de <- fread(de.file, data.table = FALSE) %>%
    mutate(gene_id = str_remove(str_trim(gene_id), "\\.t\\d+$")) %>%
    distinct(gene_id, .keep_all = TRUE)
  result <- genes %>%
    left_join(de, by = "gene_id") %>%
    mutate(direction = case_when(
      !is.na(padj) & padj < de.fdr & abs(log2FoldChange) >= de.lfc & log2FoldChange > 0 ~ "SW up",
      !is.na(padj) & padj < de.fdr & abs(log2FoldChange) >= de.lfc & log2FoldChange < 0 ~ "LW up",
      TRUE ~ "Not significant"
    ))
  message(assembly, " GTF genes matched to DESeq2: ",
          sum(!is.na(result$log2FoldChange)))
  result
}
genes <- list(
  LW = read.genes.and.de(gtf.files[["LW"]], de.files[["LW"]], "LW"),
  SW = read.genes.and.de(gtf.files[["SW"]], de.files[["SW"]], "SW")
)
gc()

# ============================================================
# Plot: mirrored LW and SW tracks around LW/SW AGP synteny
# ============================================================

all.plots <- setNames(vector("list", length(chr.order)), chr.order)
for (lw.chr in chr.order) {
  sw.chr <- as.character(best.pairs$q.name[match(lw.chr, best.pairs$t.name)])
  if (length(sw.chr) != 1 || is.na(sw.chr)) { warning("Skipping ", lw.chr); next }
  p <- paf %>% filter(t.name == lw.chr, q.name == sw.chr)
  if (!nrow(p)) { warning("No PAF data for ", lw.chr); next }
  lw.length <- max(as.numeric(p$t.len), na.rm = TRUE)
  sw.length <- max(as.numeric(p$q.len), na.rm = TRUE)
  canvas <- max(lw.length, sw.length)
  
  p.fst.lw <- make.fst.plot(
    filter(fst$LW, chrom == lw.chr, between(position, 0, lw.length)),
    canvas, "LW"
  )
  p.repeat.lw <- make.repeat.plot(
    filter(repeats$LW, seqnames == lw.chr,
           between(position, 0, lw.length)), canvas, "LW"
  )
  p.map.lw <- make.genetic.map.plot(
    filter(genetic.maps, reference_genome == "LW", chrom == lw.chr,
           between(map.position, 0, lw.length)),
    canvas, "LW"
  )
  p.deg.lw <- make.deg.plot(filter(genes$LW, chrom == lw.chr), canvas)
  p.genes.lw <- make.gene.plot(filter(genes$LW, chrom == lw.chr), canvas, "LW")
  
  p.synteny <- make.synteny.plot(
    p, filter(agp.lw, chrom == lw.chr), filter(agp.sw, chrom == sw.chr),
    lw.chr, sw.chr, lw.length, sw.length
  )
  
  p.genes.sw <- make.gene.plot(filter(genes$SW, chrom == sw.chr), canvas, "SW")
  p.deg.sw <- make.deg.plot(filter(genes$SW, chrom == sw.chr), canvas)
  p.repeat.sw <- make.repeat.plot(
    filter(repeats$SW, seqnames == sw.chr,
           between(position, 0, sw.length)), canvas, "SW"
  )
  p.map.sw <- make.genetic.map.plot(
    filter(genetic.maps, reference_genome == "SW", chrom == sw.chr,
           between(map.position, 0, sw.length)),
    canvas, "SW"
  )
  p.fst.sw <- make.fst.plot(
    filter(fst$SW, chrom == sw.chr, between(position, 0, sw.length)),
    canvas, "SW", show.x = TRUE,
    comparison.order = rev(comparison.labels)
  )
  
  # Build independent masks in each assembly's own coordinate system.
  repeat.mask.lw <- repeats$LW %>%
    filter(seqnames == lw.chr,
           normalized_coverage > repeat.mask.threshold,
           between(position, 0, lw.length)) %>%
    dplyr::select(start, end) %>%
    merge.mask.intervals()
  
  repeat.mask.sw <- repeats$SW %>%
    filter(seqnames == sw.chr,
           normalized_coverage > repeat.mask.threshold,
           between(position, 0, sw.length)) %>%
    dplyr::select(start, end) %>%
    merge.mask.intervals()
  
  # LW masks cover the upper tracks through the LW gene track.
  # They stop before the AGP scaffolds and synteny ribbons.
  p.fst.lw <- add.repeat.mask(p.fst.lw, repeat.mask.lw)
  p.map.lw <- add.repeat.mask(p.map.lw, repeat.mask.lw)
  p.repeat.lw <- add.repeat.mask(p.repeat.lw, repeat.mask.lw)
  p.deg.lw <- add.repeat.mask(p.deg.lw, repeat.mask.lw)
  p.genes.lw <- add.repeat.mask(p.genes.lw, repeat.mask.lw)
  
  # SW masks begin at the SW gene track, below the synteny panel.
  p.genes.sw <- add.repeat.mask(p.genes.sw, repeat.mask.sw)
  p.deg.sw <- add.repeat.mask(p.deg.sw, repeat.mask.sw)
  p.repeat.sw <- add.repeat.mask(p.repeat.sw, repeat.mask.sw)
  p.map.sw <- add.repeat.mask(p.map.sw, repeat.mask.sw)
  p.fst.sw <- add.repeat.mask(p.fst.sw, repeat.mask.sw)
  
  combined <- p.fst.lw / p.map.lw / p.repeat.lw / p.deg.lw / p.genes.lw /
    p.synteny / p.genes.sw / p.deg.sw / p.repeat.sw / p.map.sw / p.fst.sw +
    plot_layout(
      heights = c(2.2, 0.48, 0.24, 0.72, 0.12, 1.25,
                  0.12, 0.72, 0.24, 0.48, 2.2),
      guides = "collect"
    ) +
    plot_annotation(
      title = paste0("LW ", lw.chr, " versus SW ", sw.chr),
      subtitle = "Mirrored LW and SW FST, repeat-content and expression tracks",
      theme = theme(plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 10))
    )
  all.plots[[lw.chr]] <- combined
  prefix <- paste0("LW_", lw.chr, "_vs_SW_", sw.chr,
                   "_mirrored_FST_repeats_genes_AGP_synteny")
  ggsave(file.path(output.dir, paste0(prefix, ".pdf")), combined,
         width = 14, height = 19.5, units = "in")
  ggsave(file.path(output.dir, paste0(prefix, ".png")), combined,
         width = 14, height = 19.5, units = "in", dpi = 300)
  message("Saved ", prefix)
}

saveRDS(all.plots, file.path(output.dir, "all_LW_SW_layered_plots.rds"))
message("Finished. Output: ", output.dir)
