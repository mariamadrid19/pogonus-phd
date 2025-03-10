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
plink --bfile gwas_input --pheno phenotype.txt --assoc --out check_pheno


source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate gemma

# Generate a kinship matrix 
gemma -bfile gwas_input -gk 1 -o kinship

# Perform the GWAS with phenotype file and kinship matrix
gemma -bfile gwas_input -k output/kinship.cXX.txt -lmm 4 -o gwas_results
