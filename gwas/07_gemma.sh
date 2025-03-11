#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name gwas 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o gwas.%j.out

module load Python/3.7.0-foss-2018a
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

# Remove the samples that don't have a wing measurement (not found in the phenotype.txt file) 
bcftools view --threads 20 --samples-file ^remove_samples.txt gwas_imputed.vcf.gz -Oz -o gwas_imputed_clean.vcf.gz

bcftools index -t gwas_imputed_clean.vcf.gz

bcftools query -l gwas_imputed_clean.vcf.gz | wc -l
# total of 237 samples (some were removed due to not having a phenotype measurement)

module load PLINK/1.9
plink --vcf gwas_imputed_clean.vcf.gz --pheno phenotype_final.txt --allow-no-sex --pheno-name wingsize --double-id --make-bed --allow-extra-chr --out gwas_input
# fixed the phenotype file so that the FID and IID columns are the same, and that it is in the same order as the samples in the vcf file 

# To confirm that the .bed file is properly formatted
plink --bfile gwas_input --freq --allow-extra-chr --allow-no-sex

# This generates missing_check.imiss, which tells you if any individuals have missing genotype data
plink --bfile gwas_input --missing --out missing_check --allow-extra-chr --allow-no-sex

# Generates the genetic relationship matrix (GRM), also known as kinship matrix, which GEMMA requires -> gwas_input.grm.bin
plink --bfile gwas_input --make-grm-bin --out gwas_input --allow-extra-chr --allow-no-sex

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate gemma

# Generate a kinship matrix
gemma -bfile gwas_input -gk 1 -outdir kinship_matrix -o gwas_input

# Perform the GWAS with phenotype file and kinship matrix
gemma -bfile gwas_input -k kinship_matrix/gwas_input.cXX.txt -lmm 4 -outdir gemma_results -o gemma_lmm_results
