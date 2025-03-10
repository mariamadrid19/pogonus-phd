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

vcftools --gzvcf merged_variants.vcf.gz \
    --max-missing 0.9 \        # At least 90% of individuals must have a called genotype
    --minQ 20 \                # Minimum quality score of 30
    --maf 0.05 \               # Minor allele frequency (MAF) threshold of 0.05
    --remove-indels \          # Exclude indels, keeping only SNPs
    --recode --stdout | bgzip > gwas_filtered.vcf.gz

# Index VCF file
tabix -p vcf gwas_filtered.vcf.gz
