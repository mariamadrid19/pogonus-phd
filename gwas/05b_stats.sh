#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name stats 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o stats.%j.out

module load Python/3.7.0-foss-2018a
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

# Gives variant count, transition/transversion ratio, missing genotype rates, etc..
bcftools stats gwas_filtered.vcf.gz

# Gives a brief summary of stats
bcftools stats gwas_filtered.vcf.gz | grep -E "SN|TSTV"

# Counts the number of variants (excluding the header).
bcftools view -H gwas_filtered.vcf.gz | wc -l

# Checks missingness
vcftools --gzvcf gwas_filtered.vcf.gz --missing-indv

# Checks allele frequency distribution 
vcftools --gzvcf gwas_filtered.vcf.gz --freq
