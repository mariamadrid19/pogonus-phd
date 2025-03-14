#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name clean
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=16
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o clean.%j.out

module load Python/3.7.0-foss-2018a
module load tabix/0.2.6-GCCcore-6.4.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins
VCF=file_name

# This is not 100% necessary to do, but it does avoid a lot of headaches downstream
# Sample names should be as simple as possible, with any underscores ("_") or complicated suffixes
# Should EXACTLY match the phenotype data, same sample names and in the same order

# Extract the header, modify sample names by removing suffix and underscores
bcftools view -h $VCF.vcf.gz | awk '
BEGIN {OFS="\t"}
/^#CHROM/ {
  for (i=10; i<=NF; i++) {
    sub(".primDud.filtered.sorted.nd.bam", "", $i);  # Remove suffix (from mapping)
    gsub("_", "", $i);  # Remove underscores in sample names
  }
}
{print}
' > new_header.txt

# Reheader the VCF file with the cleaned sample names
bcftools reheader -h new_header.txt -o clean_$VCF.vcf.gz $VCF.vcf.gz
