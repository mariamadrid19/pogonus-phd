library(data.table)
library(ggplot2)

# Run this script in the directory containing the four PAF files, or change
# data_dir to that directory.
data_dir <- "/Users/mariamadrid/Downloads"

paf_files <- file.path(data_dir, c(
  "g1_g2.paf.gz",
  "g2_g3.paf.gz",
  "g3_g4.paf.gz",
  "g4_g5.paf.gz"
))

# AGP files must be in the same top-to-bottom order as genome_names.
agp_files <- file.path(data_dir, c(
  "Pchalceus_SW_CHR.fasta.agp",
  "Pchalceus_LW_CHR.fasta.agp",
  "P_olivaceus_CHR.fa.agp",
  "P_littoralis_CHR.fasta.agp",
  "P_gilvipes_CHR.fasta.agp"
))

genome_names <- c(
  "P. chalceus (SW)",
  "P. chalceus (LW)",
  "P. olivaceus",
  "P. littoralis",
  "P. gilvipes"
)

min_alignment_length <- 50000
min_mapq <- 20
# Horizontal space between consecutive chromosomes, expressed as a fraction
# of the longest chromosome in the five assemblies.
gap_fraction <- 0.05

# One track colour per assembly, from top to bottom.
assembly_colours <- c(
  "P. chalceus (SW)" = "#2E86AB",
  "P. chalceus (LW)" = "#36B37E",
  "P. olivaceus"     = "#F2B134",
  "P. littoralis"    = "#8E6CBB",
  "P. gilvipes"      = "#E76F51"
)

# Make the second shade by mixing each existing assembly colour with white.
lighten_colour <- function(colour, amount = 0.48) {
  rgb_values <- col2rgb(colour)
  rgb(
    red   = rgb_values[1L, ] + (255 - rgb_values[1L, ]) * amount,
    green = rgb_values[2L, ] + (255 - rgb_values[2L, ]) * amount,
    blue  = rgb_values[3L, ] + (255 - rgb_values[3L, ]) * amount,
    maxColorValue = 255
  )
}

assembly_shades <- data.table(
  genome_name = rep(names(assembly_colours), each = 2L),
  shade = rep(c("dark", "light"), times = length(assembly_colours)),
  track_colour = as.vector(rbind(
    unname(assembly_colours),
    lighten_colour(unname(assembly_colours))
  ))
)

# More points give smoother ribbons but larger output files.
ribbon_curve_points <- 30L

# Visual-only chromosome reversals. These transformations are applied after
# best-pair selection and do not modify the PAF files or biological results.
# Genome 1 is P. chalceus (SW).
flip_for_display <- data.table(
  genome = c(1L, 1L, 1L),
  chromosome = c("CHR3", "CHR6", "CHR11")
)

output_pdf <- file.path(data_dir, "five_genome_best_synteny_aligned.pdf")
output_png <- file.path(data_dir, "five_genome_best_synteny_aligned.png")

paf_columns <- c(
  "qname", "qlen", "qstart", "qend", "strand",
  "tname", "tlen", "tstart", "tend", "nmatch", "alnlen", "mapq"
)

read_paf <- function(filename, comparison) {
  x <- fread(filename, header = FALSE, fill = TRUE)
  if (ncol(x) < 12L) stop("Invalid PAF file: ", filename)
  setnames(x, 1:12, paf_columns)
  
  # Keep minimap2 primary records when the tp tag is present. If a PAF lacks
  # this optional tag, retain its records and rely on MAPQ/length filtering.
  if (ncol(x) > 12L) {
    tags <- as.matrix(x[, 13:ncol(x), with = FALSE])
    has_tp <- rowSums(tags == "tp:A:P", na.rm = TRUE) > 0L |
      rowSums(tags == "tp:A:S", na.rm = TRUE) > 0L
    primary <- rowSums(tags == "tp:A:P", na.rm = TRUE) > 0L
    x <- x[!has_tp | primary]
  }
  
  x[, comparison := comparison]
  x
}

paf_all <- Map(read_paf, paf_files, seq_along(paf_files))

# IMPORTANT: minimap2 target.fa query.fa writes query in PAF columns 1-4 and
# target in columns 6-9. Your commands therefore place genome i as the target
# (upper row) and genome i+1 as the query (lower row).

# Obtain chromosome lengths before alignment filtering, so a chromosome is not
# omitted merely because none of its blocks passes a plotting threshold.
length_parts <- vector("list", 8L)
k <- 1L
for (i in seq_along(paf_all)) {
  x <- paf_all[[i]]
  length_parts[[k]] <- unique(x[grepl("^CHR", tname), .(
    genome = i, chromosome = tname, length = tlen
  )])
  k <- k + 1L
  length_parts[[k]] <- unique(x[grepl("^CHR", qname), .(
    genome = i + 1L, chromosome = qname, length = qlen
  )])
  k <- k + 1L
}

chromosomes <- rbindlist(length_parts)[
  , .(length = max(length)), by = .(genome, chromosome)
]

# Natural order: CHR1, CHR2, ..., CHR11. This works for both CHR1 and
# CHR1_RagTag naming conventions.
chromosomes[, chromosome_number := suppressWarnings(as.integer(
  sub("^[^0-9]*([0-9]+).*", "\\1", chromosome)
))]
chromosomes[is.na(chromosome_number), chromosome_number := .Machine$integer.max]
setorder(chromosomes, genome, chromosome_number, chromosome)

counts <- chromosomes[, .N, by = genome]
print(counts)
if (any(counts$N != 11L)) {
  warning("Not every genome has exactly 11 CHR records; inspect the table above.")
}

chromosome_gap <- max(chromosomes$length) * gap_fraction

# Give every chromosome number one shared horizontal slot. The slot width is
# set by the longest copy among the five assemblies, so corresponding
# chromosomes start at exactly the same x coordinate without overlapping the
# following chromosome.
chromosome_slots <- chromosomes[, .(
  slot_width = max(length)
), by = chromosome_number]
setorder(chromosome_slots, chromosome_number)
chromosome_slots[, offset := cumsum(
  shift(slot_width + chromosome_gap, fill = 0)
)]
chromosomes[
  chromosome_slots,
  on = "chromosome_number",
  offset := i.offset
]
chromosomes[, `:=`(
  xmin = offset,
  xmax = offset + length,
  y = 6L - genome,
  genome_name = genome_names[genome]
)]

layout <- chromosomes[, .(genome, chromosome, offset, y)]

read_agp_scaffolds <- function(filename, genome_number) {
  # skip = "CHR" bypasses optional one-column AGP comment/header lines that
  # otherwise make fread infer a one-column table.
  x <- fread(filename, header = FALSE, sep = "\t", fill = TRUE, skip = "CHR")
  if (ncol(x) < 9L) stop("Invalid AGP file: ", filename)
  setnames(
    x, 1:9,
    c("chromosome", "object_beg", "object_end", "part_number",
      "component_type", "component_id", "component_beg",
      "component_end", "orientation")
  )
  # Only chromosome objects are plotted; some AGPs also contain unplaced
  # contigs/scaffolds as separate objects after the chromosome records.
  x <- x[grepl("^CHR", chromosome) & component_type == "W"]
  x[, `:=`(
    genome = genome_number,
    object_beg = as.numeric(object_beg),
    object_end = as.numeric(object_end),
    part_number = as.integer(part_number)
  )]
  x
}

agp_scaffolds <- rbindlist(Map(
  read_agp_scaffolds, agp_files, seq_along(agp_files)
))

# Map AGP objects to PAF chromosomes by exact length rather than trusting the
# object labels. This matters for the SW assembly: several CHR labels in its
# AGP describe different physical chromosomes than the same labels in the PAF.
agp_object_lengths <- agp_scaffolds[, .(
  length = max(object_end)
), by = .(genome, agp_chromosome = chromosome)]

chromosome_map <- merge(
  agp_object_lengths,
  chromosomes[, .(genome, chromosome, length)],
  by = c("genome", "length"),
  allow.cartesian = TRUE
)

ambiguous_map <- chromosome_map[, .N, by = .(genome, agp_chromosome)][N != 1L]
unmapped_agp <- agp_object_lengths[
  !chromosome_map, on = .(genome, agp_chromosome)
]
if (nrow(ambiguous_map) || nrow(unmapped_agp) ||
    nrow(chromosome_map) != nrow(agp_object_lengths)) {
  stop(
    "Could not map every AGP chromosome uniquely to a PAF chromosome by ",
    "exact length. Inspect chromosome_map and agp_object_lengths."
  )
}

message("AGP-to-PAF chromosome mapping:")
print(chromosome_map[order(genome, chromosome)])

agp_scaffolds[, agp_chromosome := chromosome]
agp_scaffolds[, chromosome := NULL]
agp_scaffolds <- merge(
  agp_scaffolds,
  chromosome_map[, .(genome, agp_chromosome, chromosome)],
  by = c("genome", "agp_chromosome")
)

scaffolds <- merge(
  agp_scaffolds,
  chromosomes[, .(genome, chromosome, offset, length, y, genome_name)],
  by = c("genome", "chromosome")
)
scaffolds[, chromosome_number := suppressWarnings(as.integer(
  sub("^[^0-9]*([0-9]+).*", "\\1", chromosome)
))]
setorder(scaffolds, genome, chromosome_number, chromosome, object_beg)
scaffolds[, scaffold_index := seq_len(.N), by = .(genome, chromosome)]
scaffolds[, shade := fifelse(scaffold_index %% 2L == 1L, "dark", "light")]

# AGP positions are 1-based inclusive; ggplot/PAF positions are 0-based
# half-open. Mirror boundaries for chromosomes flipped for display.
scaffolds[, `:=`(
  local_xmin = object_beg - 1,
  local_xmax = object_end
)]
scaffolds[
  flip_for_display,
  on = .(genome, chromosome),
  `:=`(
    local_xmin = length - object_end,
    local_xmax = length - (object_beg - 1)
  )
]
scaffolds[, `:=`(
  xmin = offset + local_xmin,
  xmax = offset + local_xmax
)]
scaffolds <- merge(
  scaffolds, assembly_shades,
  by = c("genome_name", "shade"),
  all.x = TRUE
)

# For every query chromosome, find the target chromosome receiving the greatest
# total aligned length. Retain all primary, high-quality blocks belonging to
# that best chromosome pair so chromosome-scale structure remains visible.
select_best_pairs <- function(x) {
  x <- x[
    grepl("^CHR", qname) & grepl("^CHR", tname) &
      alnlen >= min_alignment_length & mapq >= min_mapq
  ]
  if (!nrow(x)) return(x)
  
  scores <- x[, .(mapped_bases = sum(alnlen)), by = .(qname, tname)]
  setorder(scores, qname, -mapped_bases, tname)
  best <- scores[, .SD[1L], by = qname]
  
  message("Best chromosome pairs:")
  print(best)
  x[best, on = .(qname, tname), nomatch = 0L]
}

paf_best <- lapply(paf_all, select_best_pairs)

make_ribbons <- function(x, upper_genome, lower_genome) {
  # Work on a copy: coordinate reversals below are for plotting only.
  x <- copy(x)
  
  flipped_upper <- flip_for_display[genome == upper_genome, chromosome]
  flip_target <- x$tname %chin% flipped_upper
  if (any(flip_target)) {
    old_start <- x$tstart[flip_target]
    old_end <- x$tend[flip_target]
    x[flip_target, `:=`(
      tstart = tlen - old_end,
      tend = tlen - old_start,
      strand = fifelse(strand == "+", "-", "+")
    )]
  }
  
  flipped_lower <- flip_for_display[genome == lower_genome, chromosome]
  flip_query <- x$qname %chin% flipped_lower
  if (any(flip_query)) {
    old_start <- x$qstart[flip_query]
    old_end <- x$qend[flip_query]
    x[flip_query, `:=`(
      qstart = qlen - old_end,
      qend = qlen - old_start,
      strand = fifelse(strand == "+", "-", "+")
    )]
  }
  
  upper <- layout[genome == upper_genome, .(
    tname = chromosome, upper_offset = offset, upper_y = y
  )]
  lower <- layout[genome == lower_genome, .(
    qname = chromosome, lower_offset = offset, lower_y = y
  )]
  
  z <- merge(x, upper, by = "tname")
  z <- merge(z, lower, by = "qname")
  z[, alignment_id := sprintf("G%d_G%d_%d", upper_genome, lower_genome, .I)]
  
  # Upper sequence is the PAF target; lower sequence is the PAF query.
  z[, `:=`(
    upper_x1 = upper_offset + tstart,
    upper_x2 = upper_offset + tend,
    query_x1 = lower_offset + qstart,
    query_x2 = lower_offset + qend
  )]
  z[, `:=`(
    lower_x1 = fifelse(strand == "+", query_x1, query_x2),
    lower_x2 = fifelse(strand == "+", query_x2, query_x1)
  )]
  
  # Sample two cubic/smoothstep sides for every ribbon. Zero slope at both
  # ends makes ribbons flow smoothly out of and into the assembly tracks.
  curve_one_ribbon <- function(id, strand_value, ux1, ux2, uy, lx1, lx2, ly) {
    progress <- seq(0, 1, length.out = ribbon_curve_points)
    smooth <- 3 * progress^2 - 2 * progress^3
    
    data.table(
      alignment_id = id,
      strand = strand_value,
      vertex = seq_len(2L * ribbon_curve_points),
      x = c(
        ux1 + (lx1 - ux1) * smooth,
        lx2 + (ux2 - lx2) * smooth
      ),
      y = c(
        uy + (ly - uy) * progress,
        ly + (uy - ly) * progress
      )
    )
  }
  
  rbindlist(Map(
    curve_one_ribbon,
    z$alignment_id, z$strand,
    z$upper_x1, z$upper_x2, z$upper_y,
    z$lower_x1, z$lower_x2, z$lower_y
  ))
}

ribbons <- rbindlist(lapply(seq_along(paf_best), function(i) {
  make_ribbons(paf_best[[i]], i, i + 1L)
}))
ribbons[, orientation := fifelse(strand == "-", "Inversion", "Forward")]
setorder(ribbons, alignment_id, vertex)

chromosome_labels <- chromosomes[genome == 1L, .(
  chromosome,
  x = (xmin + xmax) / 2,
  y = y + 0.09
)]
genome_labels <- chromosomes[, .(
  x = min(xmin), y = unique(y), genome_name = unique(genome_name)
), by = genome]

p <- ggplot() +
  geom_polygon(
    data = ribbons,
    aes(x = x, y = y, group = alignment_id, fill = orientation),
    alpha = 0.35, colour = NA
  ) +
  # Rounded, thin assembly tracks with a subtle dark outline.
  geom_segment(
    data = chromosomes,
    aes(x = xmin, xend = xmax, y = y, yend = y),
    colour = "grey20", linewidth = 1.15, lineend = "round"
  ) +
  # Alternate the existing assembly colour with a lighter version for each
  # successive W/component row in the AGP. Gaps remain as the grey underlay.
  geom_segment(
    data = scaffolds,
    aes(
      x = xmin, xend = xmax, y = y, yend = y,
      colour = track_colour
    ),
    linewidth = 0.72, lineend = "butt"
  ) +
  geom_text(
    data = chromosome_labels,
    aes(x = x, y = y, label = chromosome),
    size = 2.2, angle = 45, hjust = 0
  ) +
  geom_text(
    data = genome_labels,
    aes(x = x, y = y, label = genome_name),
    hjust = 1.05, fontface = "bold", size = 3.5
  ) +
  scale_fill_manual(
    name = "Alignment orientation",
    values = c(Forward = "#3B82F6", Inversion = "#EF4444")
  ) +
  scale_colour_identity() +
  scale_x_continuous(
    labels = function(x) paste0(format(x / 1e6, trim = TRUE), "")
  ) +
  scale_y_continuous(breaks = 5:1, labels = genome_names) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Genomic position (Mb)", y = NULL,
    title = "Best chromosome mappings among five genomes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(), panel.border = element_blank(),
    axis.line.x = element_line(), axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(20, 20, 20, 110),
    legend.position = "none"
  )

ggsave(output_pdf, p, width = 22, height = 10, units = "in", limitsize = FALSE)
ggsave(output_png, p, width = 22, height = 10, units = "in", dpi = 800, limitsize = FALSE)
print(p)
