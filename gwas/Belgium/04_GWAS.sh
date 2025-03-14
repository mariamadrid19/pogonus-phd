#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name gwas 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o gwas.%j.out

module load Python/3.7.0-foss-2018a
module load tabix/0.2.6-GCCcore-6.4.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins
POPULATION=Belgium
PHENOTYPE=wing_length

module load PLINK/1.9
plink --vcf gwas_imputed_clean_$POPULATION.vcf.gz --pheno phenotype_final_$POPULATION.txt --allow-no-sex --pheno-name $PHENOTYPE --make-bed --allow-extra-chr --out gwas_input_$POPULATION
# fixed the phenotype file so that the FID and IID columns are the same, and that it is in the same order as the samples in the vcf file 

# This generates missing_check.imiss, which tells you if any individuals have missing genotype data
plink --bfile gwas_input_$POPULATION --missing --out missing_check --allow-extra-chr --allow-no-sex

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate gemma

# Generate a kinship matrix
gemma -bfile gwas_input_$POPULATION -gk 1 -outdir kinship_matrix -o gwas_input_$POPULATION

# Perform the GWAS with phenotype file and kinship matrix
gemma -bfile gwas_input_$POPULATION -k kinship_matrix/gwas_input_$POPULATION.cXX.txt -lmm 4 -outdir gemma_results -o gemma_lmm_results

# If adding the phenotype doesn't work, try this:
plink --vcf Belgium_filtered.vcf.gz --allow-extra-chr --allow-no-sex --make-bed --out Belgium_final

mv Belgium_final.fam Belgium_no_pheno.fam

# Add phenotypes with awk from phenotype file (phenos, just the quantitative values)
awk 'NR==FNR {phenotype[NR]=$1; next} {print $1, $2, $3, $4, $5, phenotype[FNR]}' phenos.txt Belgium_no_pheno.fam > Belgium_with_pheno.fam

# Add the family IDs if relevant
awk 'NR==FNR {fid[NR]=$1; next} {print fid[FNR], $2, $3, $4, $5, $6}' sample_fids.txt Belgium_with_pheno.fam > Belgium_with_pheno_updated.fam

# Add sex if relevant
awk 'NR==FNR {fid[NR]=$1; next} {print $1, $2, $3, $4, sex[FNR], $6}' sex.txt Belgium_with_pheno_updated.fam > Belgium_final.fam

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate gemma

# Generate a kinship matrix
gemma -bfile Belgium_final -gk 1 -outdir kinship_matrix -o Belgium_final

# Perform the GWAS with phenotype file and kinship matrix
gemma -bfile Belgium_final -k kinship_matrix/Belgium_final.cXX.txt -lmm 4 -outdir gemma_results -o gemma_Belgium
