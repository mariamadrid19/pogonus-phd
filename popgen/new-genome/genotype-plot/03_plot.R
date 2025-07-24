#install.packages("remotes")
#remotes::install_github("JimWhiting91/genotype_plot")
library("GenotypePlot") # requires a local installation of bcftools
setwd("/scratch/leuven/357/vsc35707/popgen/new-genome-july-2025/filtered-thinned")

# List of individuals
individuals <- c(
  "bams/GC129388.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129389.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129390.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129391.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129392.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129393.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129395.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129396.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129397.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129398.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129400.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129401.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129402.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129403.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129404.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129405.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129406.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129407.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129408.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129409.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129410.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129411.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129412.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129413.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129414.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129415.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129416.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129417.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129418.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129421.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129422.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129423.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129424.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129425.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129426.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129427.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129428.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129429.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129430.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129431.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129432.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129433.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129434.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129435.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129437.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC129438.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136078.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136079.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136080.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136081.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136082.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136083.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136084.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136085.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136086.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136087.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136088.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136089.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136090.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136091.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136092.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136093.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136094.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136095.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136096.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136097.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136098.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136099.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136100.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136101.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136102.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136103.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136104.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136107.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136108.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136109.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136110.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136111.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136112.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136113.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136114.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136115.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136116.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136117.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136118.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136119.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136120.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136121.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136122.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136123.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136124.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136125.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136126.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136127.Pchal_Bar_SW.filtered.sorted.dedup.bam",
  "bams/GC136128.Pchal_Bar_SW.filtered.sorted.dedup.bam")

# Define population assignment
pop_assignments <- c(
  rep("B_NIE_T", 6),
  rep("B_DUD_S", 4),
  rep("F_GUE_T", 6),
  rep("F_GUE_S", 6),
  rep("P_AVE_T", 5),
  rep("P_AVE_S", 3),
  rep("F_CAM_S", 2),
  rep("E_SEV_T", 3),
  rep("B_HEI_2000", 6),
  rep("B_HEI_2018", 5),
  rep("F_GUE_T", 6),
  rep("B_DUD_S", 6),
  rep("B_NIE_T", 6),
  rep("F_GUE_S", 6),
  rep("F_CAM_S", 3),
  rep("S_HUE_T", 2),
  rep("S_COT_S", 5),
  rep("S_HUE_T", 3),
  rep("B_HEI_2000", 6),
  rep("B_HEI_2018", 6)
)

# Combine into a data frame
pop_df <- data.frame(
  ind = individuals,
  #sample_id = sample_ids,
  pop = pop_assignments,
  stringsAsFactors = FALSE
)

# View the result
head(pop_df)

# Chromosomes to plot
chromosomes <- list(
  "2" = 53230543,
  "4" = 44639843,
  "5" = 65146079,
  "6" = 36339123,
  "8" = 35503798
)

# Loop through chromosomes
for (chr_num in names(chromosomes)) {
  chr_id <- paste0("CHR", chr_num)
  vcf_file <- paste0("Pchal_Bar_SW.chr_", chr_num, ".filtered.nomissing.thin500.vcf.gz")
  out_file <- paste0("genotype_plot_CHR", chr_num, ".png")
  
  cat("Processing", chr_id, "->", vcf_file, "\n")
  
  p <- genotype_plot(
    vcf = vcf_file,
    chr = chr_id,
    start = 1,
    end = chromosomes[[chr_num]],
    popmap = pop_df,
    cluster = TRUE,
    snp_label_size = 10000,
    colour_scheme = c("#FCD225", "#C92D59", "#300060"),
    invariant_filter = TRUE
  )
  
  png(out_file, width = 1600, height = 1000, res = 150)
  combine_genotype_plot(p)
  dev.off()
}
