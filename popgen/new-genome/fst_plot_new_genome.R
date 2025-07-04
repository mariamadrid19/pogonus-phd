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
setwd("/scratch/leuven/357/vsc35707/popgen")

# Load scaffold lengths
scafL <- read.table("P_chalceus_final_lengths.txt", header = FALSE, stringsAsFactors = FALSE)

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
setwd("/scratch/leuven/357/vsc35707/popgen/final-stats")
pattern_BE <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11|12)\\.pop0\\.stats"
pattern_FR <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11|12)\\.pop1\\.stats"
pattern_PO <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11|12)\\.pop2\\.stats"
pattern_SP <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11|12)\\.pop3\\.stats"
pattern_HEIST <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11|12)\\.pop4\\.stats"
pattern_UK_ME <- "Pchal_Bar_SW\\.chr_(1|2|3|4|5|6|7|8|9|10|11|12)\\.pop5\\.stats"

file_list_BE <- list.files(pattern = pattern_BE, full.names = TRUE)
file_list_FR <- list.files(pattern = pattern_FR, full.names = TRUE)
file_list_PO <- list.files(pattern = pattern_PO, full.names = TRUE)
file_list_SP <- list.files(pattern = pattern_SP, full.names = TRUE)
file_list_UM <- list.files(pattern = pattern_UK_ME, full.names = TRUE)
file_list_HE <- list.files(pattern = pattern_HEIST, full.names = TRUE)

stats_BE <- rbindlist(lapply(file_list_BE, fread))
stats_FR <- rbindlist(lapply(file_list_FR, fread))
stats_PO <- rbindlist(lapply(file_list_PO, fread))
stats_SP <- rbindlist(lapply(file_list_SP, fread))
stats_UM <- rbindlist(lapply(file_list_UM, fread))
stats_HE <- rbindlist(lapply(file_list_HE, fread))


comp=list(stats_BE,stats_FR,stats_PO,stats_SP,stats_UM,stats_HE)
names=c("Belgium T/S","France T/S","Portugal T/S", "Spain T/S", "UK (T) / FR (S)", "Heist 2000/2018")
col=c("black","black","black","black","black","black","black","black", "black")

layout(matrix(c(1:6), nrow=6, byrow=TRUE), height = c(1,1,1,1,1,1))
layout.show(n=6)

begin = chrom_coords$chromStarts[1]/1000000
end = max(chrom_coords$chromEnds)/1000000

top = 1
bot = 0
bot2 = 0


###
# plot Fst
# Set plot margins
par(mai = c(0.05, 0.4, 0.05, 0.4), oma = c(2, 0, 1, 0) + 0)

# Define color palette
colfuncR <- colorRampPalette(c("black", "red"))
colL <- colfuncR(100)

par(mai=c(0.05,0.4,0.05,0.4), oma=c(2,0,1,0)+0)
colfuncR <- colorRampPalette(c("black", "red"))
colL <- colfuncR(100)
for (i in 1:length(comp)) {
  comb <- merge(comp[[i]], scaf_coords, by = 'scaffold', all.x=TRUE)
  comb <- as.data.frame(comb)
  comb$chromPos <- comp[[i]]$mid + comb$chromStarts
  
  comb <- na.omit(comb)
  
  # Remove rows with "none" in the 16th column and convert the column to numeric
  comb <- comb %>% 
    filter(comb[, 16] != "none") %>% 
    mutate(across(16, as.numeric))
  
  # Ensure the numeric conversion worked and handle NA values if any
  comb[, 16][is.na(comb[, 16])] <- 0
  comb[,16][comb[,16] < 0] <- 0
  
  comb$col <- colfuncR(100)[as.integer((comb[,16]/1)*100)+1]
  
  plot(0, pch = "", xlim = c(begin,end), ylim = c(bot,top), ylab = "", yaxt = "n", lwd = 0.5, xlab = "", xaxt = "n", bty = "
n", main = "", axes = FALSE)
  rect(chrom_coords$chromStarts/1000000, rep(bot, length(chrom_coords$chromStarts)), chrom_coords$chromEnds/1000000, rep(top
, length(chrom_coords$chromStarts)), col = c("gray95", "gray90"), lwd = 0, border = c("gray95", "gray90"))
  
  par(new=TRUE)
  plot(comb$chromPos/1000000, comb[,16], type="p", pch=19, cex=0.7, col=adjustcolor(comb$col, alpha=1), xlim = c(begin,end),
 ylim=c(bot,top), axes = FALSE, bty = "n", xlab = "", ylab = "", yaxt="n", xaxt="n")
  
  axis(2, cex.axis = 1, line = -1.5)
  mtext(side = 2, text = expression(paste("Fst")), cex=0.5, line = 0.8)
  mtext(side = 4, text = names[i], cex=0.6)
}

axis(1, at=chrom_coords[,5][1:12]/1000000, labels=(1:12), lwd=0, lwd.ticks=0)
segments(x0=5, y0=0.9, x1=15, y1=0.9, lwd = 1)
text(10, 0.8, labels = "10Mb", cex = 1)



# Extract column 16 names
colname_FR <- colnames(stats_FR)[16]
colname_BE <- colnames(stats_BE)[16]
colname_SP <- colnames(stats_SP)[16]
colname_PO <- colnames(stats_PO)[16]

# Print and compare
cat("FR:", colname_FR, "\n")
cat("BE:", colname_BE, "\n")
cat("SP:", colname_SP, "\n")
cat("PO:", colname_PO, "\n")


comb <- merge(stats_PO, scaf_coords, by = "scaffold", all.x = TRUE)

# Calculate genomic position
comb$chromPos <- comb$mid + comb$chromStarts
comb <- na.omit(comb)

# Convert column 16 to numeric after filtering "none"
comb <- comb %>%
  filter(comb[, 16] != "none") %>%
  mutate(across(16, as.numeric))

# Count how many NAs were produced
num_NAs <- sum(is.na(comb[, 16]))

# Count how many values are negative
num_negatives <- sum(comb[, 16] < 0, na.rm = TRUE)

# Total number of rows before replacement
total_values <- nrow(comb)

# Replace invalid values with 0
comb[, 16][is.na(comb[, 16])] <- 0
comb[, 16][comb[, 16] < 0] <- 0

# Total values replaced
total_replaced <- num_NAs + num_negatives

# Output the results
cat("Converted to NA: ", num_NAs, "\n")
cat("Negative values: ", num_negatives, "\n")
cat("Total values set to zero: ", total_replaced, "\n")
cat("Total values (all rows): ", total_values, "\n")


# Replace these with your actual counts
pop_names <- c("Portugal", "Spain", "Belgium", "France")
total_snps <- c(36639, 37124, 37022, 37199)         # total SNPs
zeroed_snps <- c(31560, 14864, 23077, 14910)    # NA + negative SNPs set to zero

# Build dataframe
df <- data.frame(
  Population = rep(pop_names, each = 2),
  Category = rep(c("Total SNPs", "Set to zero"), times = 4),
  Count = c(36639, 31560, 37124, 14864, 37022, 23077, 37199, 14910)
)

# Plot
library(ggplot2)

ggplot(df, aes(x = Population, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Total SNPs vs SNPs set to zero",
       y = "Number of SNPs",
       fill = "") +
  theme_minimal(base_size = 14)
