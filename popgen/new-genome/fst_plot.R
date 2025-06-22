library(RColorBrewer)
library(dplyr)
library(ggplot2)

#calculate chromosome coordinates
chrom.coords <- function(scafL,chromNames,gap = 9000000) {
  chromosome = vector()
  chromLengths  = vector()
  chromStarts = vector()
  chromEnds = vector()
  chromMid = vector()
  chrom = 1
  endLast = 0
  scafCurrent <- subset(scafL, chromosome == chromNames[1])
  chromosome[chrom] <- chrom
  chromLengths[chrom] <- sum(scafCurrent$length)
  chromStarts[chrom] <- endLast + 1
  chromEnds[chrom] <- endLast + chromLengths[chrom]
  chromMid[chrom] <- endLast + chromLengths[chrom]/2
  endLast = chromEnds[chrom]
  chrom = chrom + 1
  for (i in 2:length(chromNames)) {
    chromosome[chrom] <- chrom
    scafCurrent <- subset(scafL, chromosome == chromNames[i])
    chromLengths[chrom] <- sum(scafCurrent$length)
    chromStarts[chrom] <- endLast + gap + 1
    chromEnds[chrom] <- endLast + gap + chromLengths[chrom]
    chromMid[chrom] <- endLast + gap + chromLengths[chrom]/2
    endLast = chromEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(chromosome,chromLengths,chromStarts,chromEnds,chromMid))
  return(table)
}

#calculate scaffold coordinates
scaf.coords <- function(scafL,gap = 0) {
  scaffold = vector()
  scafStarts = vector()
  scafEnds = vector()
  chrom = 1
  endLast = 0
  for (e in 1:nrow(scafL)) {
    scaffold[chrom] <- levels(scafL$scaffold)[e]
    scafStarts[chrom] <- endLast + gap + 1
    scafEnds[chrom] <- endLast + gap + scafL$length[e]
    endLast = scafEnds[chrom]
    chrom = chrom + 1
  }
  table <- as.data.frame(cbind(scaffold,scafStarts,scafEnds))
  return(table)
}

setwd("/scratch/leuven/357/vsc35707/popgen")

scafL <- read.table("new_genome_lengths.txt", h=T)
chromNames <-c(1:10)
chrom_coords <- chrom.coords(scafL, chromNames)
chrom_coords

scaf_coords <- merge(scafL, chrom_coords, by="chromosome", all.x=TRUE)
write.table(scaf_coords, "scaf_coords.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


## NEW STATS, NEW GENOME
setwd("/scratch/leuven/357/vsc35707/popgen/final-stats")
pattern_BE <- "Pchal\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.pop0\\.stats"
pattern_SP <- "Pchal\\.chr_(1|2|3|4|5|6|7|8|9|10)\\.pop1\\.stats"
file_list_BE <- list.files(pattern = pattern_BE, full.names = TRUE)
file_list_SP <- list.files(pattern = pattern_SP, full.names = TRUE)
stats_BE <- rbindlist(lapply(file_list_BE, fread))
stats_SP <- rbindlist(lapply(file_list_SP, fread))

## PLOT SPANISH STATS
comb <- merge(stats_SP, scaf_coords, by = "scaffold", all.x = TRUE)
comb$chromPos <- comb$mid + comb$chromStarts
comb <- na.omit(comb)

# Clean Fst values
comb <- comb %>%
  filter(Fst_S_HUE_T_S_COT_S != "none") %>%
  mutate(Fst_S_HUE_T_S_COT_S = as.numeric(Fst_S_HUE_T_S_COT_S))

comb$Fst_S_HUE_T_S_COT_S[is.na(comb$Fst_S_HUE_T_S_COT_S)] <- 0
comb$Fst_S_HUE_T_S_COT_S[comb$Fst_S_HUE_T_S_COT_S < 0] <- 0

# Color scale
colfunc <- colorRampPalette(c("black", "red"))
comb$col <- colfunc(100)[as.integer((comb$Fst_S_HUE_T_S_COT_S / 0.8) * 100) + 1]

# Plot settings
begin <- chrom_coords$chromStarts[1] / 1e6
end   <- chrom_coords$chromEnds[10] / 1e6
top   <- 1
bot   <- 0

# Plot
par(mfrow = c(1,1), mai = c(1,1,0.5,0.5))
plot(0, type = "n", xlim = c(begin, end), ylim = c(bot, top),
     xlab = "Chromosome", ylab = expression(F[ST]), xaxt = "n", yaxt = "s", bty = "n")

# Background alternating chromosome blocks
rect(chrom_coords$chromStarts / 1e6, rep(bot, 10),
     chrom_coords$chromEnds / 1e6, rep(top, 10),
     col = rep(c("gray95", "gray90"), length.out = 10), border = NA)

# Overlay Fst points
points(comb$chromPos / 1e6, comb$Fst_S_HUE_T_S_COT_S,
       pch = 19, cex = 0.6, col = comb$col)

# X-axis: chromosome numbers
axis(1, at = chrom_coords$chromStarts / 1e6 + (chrom_coords$chromEnds - chrom_coords$chromStarts) / 2e6,
     labels = 1:10, tick = FALSE)


# Merge with scaffold coordinate
comb <- merge(stats_BE, scaf_coords, by = "scaffold", all.x = TRUE)
comb$chromPos <- comb$mid + comb$chromStarts
comb <- na.omit(comb)

# Clean Fst values
comb <- comb %>%
  filter(Fst_B_NIE_T_B_DUD_S != "none") %>%
  mutate(Fst_B_NIE_T_B_DUD_S = as.numeric(Fst_B_NIE_T_B_DUD_S))

comb$Fst_B_NIE_T_B_DUD_S[is.na(comb$Fst_B_NIE_T_B_DUD_S)] <- 0
comb$Fst_B_NIE_T_B_DUD_S[comb$Fst_B_NIE_T_B_DUD_S < 0] <- 0

# Color scale (black to red)
colfunc <- colorRampPalette(c("black", "red"))
comb$col <- colfunc(100)[as.integer((comb$Fst_B_NIE_T_B_DUD_S / 0.8) * 100) + 1]

# Plot settings
begin <- chrom_coords$chromStarts[1] / 1e6
end   <- chrom_coords$chromEnds[10] / 1e6
top   <- 1
bot   <- 0

# Plot
par(mfrow = c(1,1), mai = c(1,1,0.5,0.5))
plot(0, type = "n", xlim = c(begin, end), ylim = c(bot, top),
     xlab = "Chromosome", ylab = expression(F[ST]), xaxt = "n", yaxt = "s", bty = "n")

# Background alternating chromosome blocks
rect(chrom_coords$chromStarts / 1e6, rep(bot, 10),
     chrom_coords$chromEnds / 1e6, rep(top, 10),
     col = rep(c("gray95", "gray90"), length.out = 10), border = NA)

# Overlay Fst points
points(comb$chromPos / 1e6, comb$Fst_B_NIE_T_B_DUD_S,
       pch = 19, cex = 0.6, col = comb$col)

# X-axis: chromosome numbers
axis(1, at = chrom_coords$chromStarts / 1e6 + (chrom_coords$chromEnds - chrom_coords$chromStarts) / 2e6,
     labels = 1:10, tick = FALSE)

## PLOT BELGIUM STATS
# Merge with scaffold coordinates (assume scaf_coords and chrom_coords already loaded)
comb <- merge(stats_BE, scaf_coords, by = "scaffold", all.x = TRUE)
comb$chromPos <- comb$mid + comb$chromStarts
comb <- na.omit(comb)

# Clean Fst values
comb <- comb %>%
  filter(Fst_B_NIE_T_B_DUD_S != "none") %>%
  mutate(Fst_B_NIE_T_B_DUD_S = as.numeric(Fst_B_NIE_T_B_DUD_S))

comb$Fst_B_NIE_T_B_DUD_S[is.na(comb$Fst_B_NIE_T_B_DUD_S)] <- 0
comb$Fst_B_NIE_T_B_DUD_S[comb$Fst_B_NIE_T_B_DUD_S < 0] <- 0

# Color scale (black → red)
colfunc <- colorRampPalette(c("black", "red"))
comb$col <- colfunc(100)[as.integer((comb$Fst_B_NIE_T_B_DUD_S / 0.8) * 100) + 1]

# Plot settings
begin <- chrom_coords$chromStarts[1] / 1e6
end   <- chrom_coords$chromEnds[10] / 1e6
top   <- 1
bot   <- 0

# Plot
par(mfrow = c(1,1), mai = c(1,1,0.5,0.5))
plot(0, type = "n", xlim = c(begin, end), ylim = c(bot, top),
     xlab = "Chromosome", ylab = expression(F[ST]), xaxt = "n", yaxt = "s", bty = "n")

# Background alternating chromosome blocks
rect(chrom_coords$chromStarts / 1e6, rep(bot, 10),
     chrom_coords$chromEnds / 1e6, rep(top, 10),
     col = rep(c("gray95", "gray90"), length.out = 10), border = NA)

# Overlay Fst points
points(comb$chromPos / 1e6, comb$Fst_B_NIE_T_B_DUD_S,
       pch = 19, cex = 0.6, col = comb$col)

# X-axis: chromosome numbers
axis(1, at = chrom_coords$chromStarts / 1e6 + (chrom_coords$chromEnds - chrom_coords$chromStarts) / 2e6,
     labels = 1:10, tick = FALSE)





### PLOT BOTH
# Merge scaffold coordinates
comb_SP <- merge(stats_SP, scaf_coords, by = "scaffold", all.x = TRUE)
comb_SP$chromPos <- comb_SP$mid + comb_SP$chromStarts
comb_SP <- na.omit(comb_SP)

comb_BE <- merge(stats_BE, scaf_coords, by = "scaffold", all.x = TRUE)
comb_BE$chromPos <- comb_BE$mid + comb_BE$chromStarts
comb_BE <- na.omit(comb_BE)

# Clean Fst values
comb_SP <- comb_SP %>%
  filter(Fst_S_HUE_T_S_COT_S != "none") %>%
  mutate(Fst_S_HUE_T_S_COT_S = as.numeric(Fst_S_HUE_T_S_COT_S))
comb_SP$Fst_S_HUE_T_S_COT_S[is.na(comb_SP$Fst_S_HUE_T_S_COT_S)] <- 0
comb_SP$Fst_S_HUE_T_S_COT_S[comb_SP$Fst_S_HUE_T_S_COT_S < 0] <- 0

comb_BE <- comb_BE %>%
  filter(Fst_B_NIE_T_B_DUD_S != "none") %>%
  mutate(Fst_B_NIE_T_B_DUD_S = as.numeric(Fst_B_NIE_T_B_DUD_S))
comb_BE$Fst_B_NIE_T_B_DUD_S[is.na(comb_BE$Fst_B_NIE_T_B_DUD_S)] <- 0
comb_BE$Fst_B_NIE_T_B_DUD_S[comb_BE$Fst_B_NIE_T_B_DUD_S < 0] <- 0

# Colors
colfunc <- colorRampPalette(c("black", "red"))
comb_SP$col <- colfunc(100)[as.integer((comb_SP$Fst_S_HUE_T_S_COT_S / 0.8) * 100) + 1]
comb_BE$col <- colfunc(100)[as.integer((comb_BE$Fst_B_NIE_T_B_DUD_S / 0.8) * 100) + 1]

# Plot settings
begin <- chrom_coords$chromStarts[1] / 1e6
end   <- chrom_coords$chromEnds[10] / 1e6
top   <- 1
bot   <- 0

midpoints <- (chrom_coords$chromStarts + chrom_coords$chromEnds) / 2 / 1e6

# ======================
# Output as PNG
png("new_genome_fst.png", width = 2000, height = 1200, res = 300)
par(mfrow = c(2, 1), mar = c(2, 4.5, 1, 4), oma = c(3, 0, 1, 0))

# Top plot: Spain
plot(0, type = "n", xlim = c(begin, end), ylim = c(bot, top),
     xlab = "", ylab = expression(F[ST]), xaxt = "n", yaxt = "s", bty = "n")
rect(chrom_coords$chromStarts / 1e6, rep(bot, 10),
     chrom_coords$chromEnds / 1e6, rep(top, 10),
     col = rep(c("gray95", "gray90"), length.out = 10), border = NA)
points(comb_SP$chromPos / 1e6, comb_SP$Fst_S_HUE_T_S_COT_S,
       pch = 19, cex = 0.6, col = comb_SP$col)
mtext("Spain (SW/LW)", side = 4, line = 2.5, cex = 0.9)

# Bottom plot: Belgium
plot(0, type = "n", xlim = c(begin, end), ylim = c(bot, top),
     xlab = "", ylab = expression(F[ST]), xaxt = "n", yaxt = "s", bty = "n")
rect(chrom_coords$chromStarts / 1e6, rep(bot, 10),
     chrom_coords$chromEnds / 1e6, rep(top, 10),
     col = rep(c("gray95", "gray90"), length.out = 10), border = NA)
points(comb_BE$chromPos / 1e6, comb_BE$Fst_B_NIE_T_B_DUD_S,
       pch = 19, cex = 0.6, col = comb_BE$col)
mtext("Belgium (SW/LW)", side = 4, line = 2.5, cex = 0.9)

# Shared X-axis
axis(1, at = midpoints, labels = 1:10, tick = FALSE, line = 1)
mtext("Chromosome", side = 1, outer = TRUE, line = 1.8, cex = 1)

dev.off()

# ======================
# Output as PDF
pdf("new_genome_fst.pdf", width = 8, height = 6)
par(mfrow = c(2, 1), mar = c(2, 4.5, 1, 4), oma = c(3, 0, 1, 0))

# Top plot: Spain
plot(0, type = "n", xlim = c(begin, end), ylim = c(bot, top),
     xlab = "", ylab = expression(F[ST]), xaxt = "n", yaxt = "s", bty = "n")
rect(chrom_coords$chromStarts / 1e6, rep(bot, 10),
     chrom_coords$chromEnds / 1e6, rep(top, 10),
     col = rep(c("gray95", "gray90"), length.out = 10), border = NA)
points(comb_SP$chromPos / 1e6, comb_SP$Fst_S_HUE_T_S_COT_S,
       pch = 19, cex = 0.6, col = comb_SP$col)
mtext("Spain (SW/LW)", side = 4, line = 2.5, cex = 0.9)

# Bottom plot: Belgium
plot(0, type = "n", xlim = c(begin, end), ylim = c(bot, top),
     xlab = "", ylab = expression(F[ST]), xaxt = "n", yaxt = "s", bty = "n")
rect(chrom_coords$chromStarts / 1e6, rep(bot, 10),
     chrom_coords$chromEnds / 1e6, rep(top, 10),
     col = rep(c("gray95", "gray90"), length.out = 10), border = NA)
points(comb_BE$chromPos / 1e6, comb_BE$Fst_B_NIE_T_B_DUD_S,
       pch = 19, cex = 0.6, col = comb_BE$col)
mtext("Belgium (SW/LW)", side = 4, line = 2.5, cex = 0.9)

# Shared X-axis
axis(1, at = midpoints, labels = 1:10, tick = FALSE, line = 1)
mtext("Chromosome", side = 1, outer = TRUE, line = 1.8, cex = 1)

dev.off()



# SAVE THE DATA
# Save Spain Fst data
data.table::fwrite(comb_SP, "new_genome_stats_Spain.tsv", sep = "\t")

# Save Belgium Fst data
data.table::fwrite(comb_BE, "new_genome_stats_Belgium.tsv", sep = "\t")

data.table::fwrite(chrom_coords, "chromosome_coordinates.tsv", sep = "\t")
