#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name check_names 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o check_names.%j.out

module load Python/3.7.0-foss-2018a
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

bcftools query -l merged_variants.vcf.gz > vcf_samples.txt
cut -f2 phenotype_modified.txt  > phenotype_samples.txt
comm -3 <(sort vcf_samples.txt) <(sort phenotype_samples.txt)
