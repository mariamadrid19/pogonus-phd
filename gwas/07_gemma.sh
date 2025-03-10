#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name gwas 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o gwas.%j.out

module load PLINK/1.9
plink --vcf gwas_imputed.vcf.gz --pheno phenotype.txt --make-bed --out gwas_input

# confirm that PLINK correctly read the phenotype file
plink --bfile gwas_imputed_with_pheno --pheno phenotype.txt --assoc --out check_pheno
