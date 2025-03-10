#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name filter 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o filter_vcf.%j.out

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate vcftools
module load tabix/0.2.6-GCCcore-6.4.0

# Check missingness per individual
vcftools --gzvcf merged_variants.vcf.gz --missing-indv

# Filter out samples with more than 10% missing genotypes (keep those with ≥90% call rate)
awk '$5 <= 0.1 {print $1}' out.imiss > high_call_samples.txt

# Use --keep to retain only these samples
vcftools --gzvcf merged_variants.vcf.gz --keep high_call_samples.txt --recode --stdout | bgzip > filtered_merged_variants.vcf.gz

# filter sites by quality scores
vcftools --gzvcf filtered_merged_variants.vcf.gz --minGQ 20 --recode --stdout | bgzip > quality_filtered_merged_variants.vcf.gz

# remove indels, keep only snps
vcftools --gzvcf quality_filtered_merged_variants.vcf.gz --remove-indels --recode --stdout | bgzip > snps_only.vcf.gz

# MAF filtering
vcftools --gzvcf snps_only.vcf.gz --maf 0.05 --recode --stdout | bgzip > gwas_filtered.vcf.gz

# Index VCF file
tabix -p vcf gwas_filtered.vcf.gz
