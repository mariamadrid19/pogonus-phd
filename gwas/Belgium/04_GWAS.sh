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
POPULATION=Belgium
PHENOTYPE=wing_length

# Change sample names, avoid underscores in the sample names to avoid problems with GEMMA downstream
bcftools view gwas_imputed_Belgium.vcf.gz | \
awk '{if(NR==1 && substr($1,1,1)=="#") {for (i=10; i<=NF; i++) gsub("_", "", $i)} print}' | \
bgzip > gwas_imputed_clean_$POPULATION.vcf.gz
tabix -p vcf gwas_imputed_clean_$POPULATION.vcf.gz

bcftools query -l gwas_imputed_clean_$POPULATION.vcf.gz | wc -l

module load PLINK/1.9
plink --vcf gwas_imputed_clean_$POPULATION.vcf.gz --pheno phenotype_final_$POPULATION.txt --allow-no-sex --pheno-name $PHENOTYPE --double-id --make-bed --allow-extra-chr --out gwas_input_$POPULATION
# fixed the phenotype file so that the FID and IID columns are the same, and that it is in the same order as the samples in the vcf file 

# To confirm that the .bed file is properly formatted
plink --bfile gwas_input_$POPULATION --freq --allow-extra-chr --allow-no-sex

# This generates missing_check.imiss, which tells you if any individuals have missing genotype data
plink --bfile gwas_input_$POPULATION --missing --out missing_check --allow-extra-chr --allow-no-sex

# Generates the genetic relationship matrix (GRM), also known as kinship matrix, which GEMMA requires -> gwas_input.grm.bin
plink --bfile gwas_input_$POPULATION --make-grm-bin --out gwas_input_$POPULATION --allow-extra-chr --allow-no-sex

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate gemma

# Generate a kinship matrix
gemma -bfile gwas_input_$POPULATION -gk 1 -outdir kinship_matrix -o gwas_input_$POPULATION

# Perform the GWAS with phenotype file and kinship matrix
gemma -bfile gwas_input_$POPULATION -k kinship_matrix/gwas_input_$POPULATION.cXX.txt -lmm 4 -outdir gemma_results -o gemma_lmm_results
