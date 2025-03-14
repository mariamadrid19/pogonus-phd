#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name merge_BE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o merge_BE.%j.out

module load Python/3.7.0-foss-2018a
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

POPULATION=Belgium

cd /scratch/leuven/357/vsc35707/GWAS/$POPULATION

bcftools concat -Oz -o merged_$POPULATION.vcf.gz *.vcf.gz

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate vcftools
module load tabix/0.2.6-GCCcore-6.4.0

# Check missingness per individual
vcftools --gzvcf merged_$POPULATION.vcf.gz --missing-indv

vcftools --gzvcf merged_$POPULATION.vcf.gz --max-missing 0.9 --minQ 30 --maf 0.05 --remove-indels --recode --stdout | bgzip > gwas_filtered_$POPULATION.vcf.gz
# At least 90% of individuals must have a called genotype
# Minimum quality score of 30
# Minor allele frequency (MAF) threshold of 0.05
# Exclude indels, keeping only SNPs

# Index VCF file
bcftools index -t vcf gwas_filtered_$POPULATION.vcf.gz
