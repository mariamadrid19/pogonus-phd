library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)

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

# Set working directory
setwd("/scratch/leuven/357/vsc35707/popgen/new-genome-july-2025")

# Load scaffold lengths
scafL <- read.table("Pchalceus_SW.sorted.lengths.txt", header = FALSE, stringsAsFactors = FALSE)

# Add column names
colnames(scafL) <- c("chrom_name", "length")

# Extract CHR entries only
scafL <- scafL[grepl("^CHR", scafL$chrom_name), ]

# Make numeric
scafL$chrom_id <- as.numeric(gsub("CHR", "", scafL$chrom_name))

# Prepare input for function
scafL$chromosome <- scafL$chrom_id

# Make numeric
scafL$length <- as.numeric(scafL$length)

# Run function
chrom_coords <- chrom.coords(scafL = scafL, chromNames = 1:12)
print(chrom_coords)

# Merge chromosome and scaffold datasets
scaf_coords <- merge(scafL, chrom_coords, by="chromosome", all.x=TRUE)
names(scaf_coords)[names(scaf_coords) == 'chrom_name'] <- 'scaffold'
scaf_coords$scaffold<-as.character(scaf_coords$scaffold)


### ============ NEW GENOME, NEW STATS, NEW PLOTS 01-07-2025 =============
setwd("/scratch/leuven/357/vsc35707/popgen/new-genome-july-2025/final-stats")
pattern_BE <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop0\\.stats"
pattern_FR <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop1\\.stats"
pattern_PO <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop2\\.stats"
pattern_SP <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop3\\.stats"
pattern_UK_ME <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop4\\.stats"

file_list_BE <- list.files(pattern = pattern_BE, full.names = TRUE)
file_list_FR <- list.files(pattern = pattern_FR, full.names = TRUE)
file_list_PO <- list.files(pattern = pattern_PO, full.names = TRUE)
file_list_SP <- list.files(pattern = pattern_SP, full.names = TRUE)
file_list_UM <- list.files(pattern = pattern_UK_ME, full.names = TRUE)

stats_BE <- rbindlist(lapply(file_list_BE, fread))
stats_FR <- rbindlist(lapply(file_list_FR, fread))
stats_PO <- rbindlist(lapply(file_list_PO, fread))
stats_SP <- rbindlist(lapply(file_list_SP, fread))
stats_UM <- rbindlist(lapply(file_list_UM, fread))


comp=list(stats_BE,stats_FR,stats_PO,stats_SP,stats_UM)
names=c("Belgium T/S","France T/S","Portugal T/S", "Spain T/S", "UK (T) / FR (S)")
col=c("black","black","black","black","black","black","black","black")

layout(matrix(c(1:5), nrow=5, byrow=TRUE), height = c(1,1,1,1,1))
layout.show(n=5)

begin = chrom_coords$chromStarts[1]/1000000
end = max(chrom_coords$chromEnds)/1000000

# Set top and bottom plot limits
top <- 1
bot <- 0

# Number of chromosomes
num_chr <- 11

# Compute chromosome midpoints for labeling
chrom_coords$chromMids <- (chrom_coords$chromStarts + chrom_coords$chromEnds) / 2

# Color function for points
colfuncR <- colorRampPalette(c("black", "red"))
colL <- colfuncR(100)

# Set outer margins
par(mai = c(0.05, 0.4, 0.05, 0.4), oma = c(2, 0, 1, 0))

# Loop over all comparisons
for (i in 1:length(comp)) {
  # Merge comparison data with chromosome coordinates
  comb <- merge(comp[[i]], scaf_coords, by = 'scaffold', all.x = TRUE)
  comb <- as.data.frame(comb)
  comb$chromPos <- comb$mid + comb$chromStarts
  comb <- na.omit(comb)
  
  # Clean up the Fst column (assumed to be the 16th)
  comb <- comb %>% 
    filter(comb[, 16] != "none") %>%
    mutate(across(16, as.numeric))
  
  comb[, 16][is.na(comb[, 16])] <- 0
  comb[, 16][comb[, 16] < 0] <- 0
  comb$col <- colfuncR(100)[as.integer((comb[, 16] / 1) * 100) + 1]
  
  # Set up an empty plot
  plot(0, pch = "", xlim = c(min(chrom_coords$chromStarts), max(chrom_coords$chromEnds)) / 1e6,
       ylim = c(bot, top), ylab = "", yaxt = "n", lwd = 0.5,
       xlab = "", xaxt = "n", bty = "n", main = "", axes = FALSE)
  
  # Add alternating gray background per chromosome
  rect(chrom_coords$chromStarts / 1e6, rep(bot, num_chr),
       chrom_coords$chromEnds / 1e6, rep(top, num_chr),
       col = rep(c("gray95", "gray90"), length.out = num_chr),
       border = NA)
  
  # Overlay Fst points
  par(new = TRUE)
  plot(comb$chromPos / 1e6, comb[, 16], type = "p", pch = 19, cex = 0.7,
       col = adjustcolor(comb$col, alpha.f = 1),
       xlim = c(min(chrom_coords$chromStarts), max(chrom_coords$chromEnds)) / 1e6,
       ylim = c(bot, top), axes = FALSE, bty = "n", xlab = "", ylab = "")
  
  # Add y-axis
  axis(2, cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(F[ST]), cex = 0.5, line = 0.8)
  
  # Add comparison label
  mtext(side = 4, text = names[i], cex = 0.6)
}

# Add chromosome number labels on x-axis
axis(1, at = chrom_coords$chromMids[1:num_chr] / 1e6,
     labels = 1:num_chr, lwd = 0, lwd.ticks = 0)

# Add scale bar
segments(x0 = 5, y0 = 0.9, x1 = 15, y1 = 0.9, lwd = 1)
text(10, 0.8, labels = "10Mb", cex = 1)


### ============ SUBSAMPLE 23-07-2025 =============
setwd("/scratch/leuven/357/vsc35707/popgen/new-genome-july-2025/subsample/final-stats")
pattern_BE <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop0\\.stats"
pattern_FR <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop1\\.stats"
pattern_PO <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop2\\.stats"
pattern_SP <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop3\\.stats"
pattern_HE <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11)\\.pop4\\.stats"

file_list_BE <- list.files(pattern = pattern_BE, full.names = TRUE)
file_list_FR <- list.files(pattern = pattern_FR, full.names = TRUE)
file_list_PO <- list.files(pattern = pattern_PO, full.names = TRUE)
file_list_SP <- list.files(pattern = pattern_SP, full.names = TRUE)
file_list_HE <- list.files(pattern = pattern_HE, full.names = TRUE)

stats_BE <- rbindlist(lapply(file_list_BE, fread))
stats_FR <- rbindlist(lapply(file_list_FR, fread))
stats_PO <- rbindlist(lapply(file_list_PO, fread))
stats_SP <- rbindlist(lapply(file_list_SP, fread))
stats_HE <- rbindlist(lapply(file_list_HE, fread))


comp=list(stats_BE,stats_FR,stats_PO,stats_SP,stats_HE)
names=c("Belgium T/S","France T/S","Portugal T/S", "Spain T/S", "Heist 2000/2018")
col=c("black","black","black","black","black","black","black","black")

layout(matrix(c(1:5), nrow=5, byrow=TRUE), height = c(1,1,1,1,1))
layout.show(n=5)

begin = chrom_coords$chromStarts[1]/1000000
end = max(chrom_coords$chromEnds)/1000000

# Set top and bottom plot limits
top <- 1
bot <- 0

# Number of chromosomes
num_chr <- 11

# Compute chromosome midpoints for labeling
chrom_coords$chromMids <- (chrom_coords$chromStarts + chrom_coords$chromEnds) / 2

# Color function for points
colfuncR <- colorRampPalette(c("black", "red"))
colL <- colfuncR(100)

# Set outer margins
par(mai = c(0.05, 0.4, 0.05, 0.4), oma = c(2, 0, 1, 0))

# Loop over all comparisons
for (i in 1:length(comp)) {
  # Merge comparison data with chromosome coordinates
  comb <- merge(comp[[i]], scaf_coords, by = 'scaffold', all.x = TRUE)
  comb <- as.data.frame(comb)
  comb$chromPos <- comb$mid + comb$chromStarts
  comb <- na.omit(comb)
  
  # Clean up the Fst column (assumed to be the 16th)
  comb <- comb %>% 
    filter(comb[, 16] != "none") %>%
    mutate(across(16, as.numeric))
  
  comb[, 16][is.na(comb[, 16])] <- 0
  comb[, 16][comb[, 16] < 0] <- 0
  comb$col <- colfuncR(100)[as.integer((comb[, 16] / 1) * 100) + 1]
  
  # Set up an empty plot
  plot(0, pch = "", xlim = c(min(chrom_coords$chromStarts), max(chrom_coords$chromEnds)) / 1e6,
       ylim = c(bot, top), ylab = "", yaxt = "n", lwd = 0.5,
       xlab = "", xaxt = "n", bty = "n", main = "", axes = FALSE)
  
  # Add alternating gray background per chromosome
  rect(chrom_coords$chromStarts / 1e6, rep(bot, num_chr),
       chrom_coords$chromEnds / 1e6, rep(top, num_chr),
       col = rep(c("gray95", "gray90"), length.out = num_chr),
       border = NA)
  
  # Overlay Fst points
  par(new = TRUE)
  plot(comb$chromPos / 1e6, comb[, 16], type = "p", pch = 19, cex = 0.7,
       col = adjustcolor(comb$col, alpha.f = 1),
       xlim = c(min(chrom_coords$chromStarts), max(chrom_coords$chromEnds)) / 1e6,
       ylim = c(bot, top), axes = FALSE, bty = "n", xlab = "", ylab = "")
  
  # Add y-axis
  axis(2, cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(F[ST]), cex = 0.5, line = 0.8)
  
  # Add comparison label
  mtext(side = 4, text = names[i], cex = 0.6)
}

# Add chromosome number labels on x-axis
axis(1, at = chrom_coords$chromMids[1:num_chr] / 1e6,
     labels = 1:num_chr, lwd = 0, lwd.ticks = 0)

# Add scale bar
segments(x0 = 5, y0 = 0.9, x1 = 15, y1 = 0.9, lwd = 1)
text(10, 0.8, labels = "10Mb", cex = 1)
