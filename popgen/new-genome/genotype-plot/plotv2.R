library(GenotypePlot)
library(vcfR)
setwd("~/Downloads")

Sys.setenv(PATH = paste("/Users/mariamadrid/bcftools", Sys.getenv("PATH"), sep=":"))

vcf_file <- "Belgium_Pchal_Bar_SW.filtered.multiSplit.renamed.vcf.gz"

# Read once just to get sample names
vcf <- read.vcfR(vcf_file)
individuals <- colnames(vcf@gt)[-1]

pop_assignments <- ifelse(
  grepl("^Pc25Np", individuals), "Nieuwpoort",
  ifelse(grepl("^Pc25Dud", individuals), "Dudzele", NA)
)

pop_df <- data.frame(
  ind = individuals,
  pop = pop_assignments,
  stringsAsFactors = FALSE
)

print(head(pop_df))
print(table(pop_df$pop, useNA = "ifany"))

chromosomes <- c(
  CHR1 = 53807401,
  CHR2 = 53230543,
  CHR3 = 46351446,
  CHR4 = 44639843,
  CHR5 = 65146079,
  CHR6 = 36339123,
  CHR7 = 47275045,
  CHR8 = 35503798,
  CHR9 = 27201578,
  CHR10 = 22848598,
  CHR11 = 49412164
)

dir.create("genotype_plots", showWarnings = FALSE)

for (chr_id in names(chromosomes)) {
  out_file <- file.path("genotype_plots", paste0("genotype_plot_", chr_id, ".png"))
  cat("Processing", chr_id, "->", out_file, "\n")
  
  p <- genotype_plot(
    vcf = vcf_file,
    chr = chr_id,
    start = 1,
    end = chromosomes[[chr_id]],
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
