#!/bin/bash -l 
#SBATCH --cluster=wice
#SBATCH --job-name gwas 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o gwas.%j.out

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh

conda activate variant_tools

# Phenotype file should have the same FID and IID columns, and sample order should be same in vcf file 
module load PLINK/1.9
plink --vcf final-vcfs/Pchal_Bar_SW.filtered.multiSplit.renamed.vcf.gz --pheno phenotypes.txt --allow-no-sex --pheno-name relMRWS --double-id --make-bed --allow-extra-chr --out gwas_input

# To confirm that the .bed file is properly formatted
plink --bfile gwas_input --freq --allow-extra-chr --allow-no-sex

# This generates missing_check.imiss, which tells you if any individuals have missing genotype data
plink --bfile gwas_input --missing --out missing_check --allow-extra-chr --allow-no-sex

conda activate gemma

# Generate a kinship matrix
gemma -bfile gwas_input -gk 1 -outdir kinship_matrix -o gwas_input

# Perform the GWAS with phenotype file and kinship matrix
# -lmm 4 linear mixed model doing p-wald, likelihood and score tests
gemma -bfile gwas_input -k kinship_matrix/gwas_input.cXX.txt -lmm 4 -outdir gemma_results_relMRWS -o gemma_lmm_results_relMRWS
