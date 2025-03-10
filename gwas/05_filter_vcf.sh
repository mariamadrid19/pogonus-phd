#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name filter 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o filter_vcf.%j.out

module load Python/3.7.0-foss-2018a
module load VCFtools/0.1.16-GCC-10.3.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

# filter sites by call rate (missing data)
vcftools --vcf merged_variants.vcf --max-missing 0.9 --recode --stdout | bgzip > filtered_merged_variants.vcf.gz

# filter sites by quality scores
vcftools --vcf filtered_merged_variants.vcf.gz --minGQ 20 --recode --stdout | bgzip > quality_filtered_merged_variants.vcf.gz

# remove indels, keep only snps
vcftools --vcf quality_filtered_merged_variants.vcf.gz --remove-indels --recode --stdout | bgzip > snps_only.vcf.gz

# MAF filtering
vcftools --vcf snps_only.vcf.gz --maf 0.05 --recode --stdout | bgzip > gwas_filtered.vcf.gz

# Index VCF file
module load tabix/0.2.6-GCCcore-6.4.0
tabix -p vcf gwas_filtered.vcf.gz
