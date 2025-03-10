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
bcftools index gwas_filtered_clean.vcf.gz

module load PLINK/1.9
plink --vcf gwas_imputed_clean.vcf.gz --pheno phenotype.txt --make-bed --out gwas_input

# confirm that PLINK correctly read the phenotype file
plink --bfile gwas_input --pheno phenotype.txt --assoc --out check_pheno

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate gemma

# Generate a kinship matrix 
gemma -bfile gwas_input -gk 1 -o kinship

# Perform the GWAS with phenotype file and kinship matrix
gemma -bfile gwas_input -k output/kinship.cXX.txt -lmm 4 -o gwas_results
