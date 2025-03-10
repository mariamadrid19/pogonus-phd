#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name impute 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o impute.%j.out

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate vcftools
module load Python/3.7.0-foss-2018a
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

# Beagle 5.5
zcat gwas_filtered.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > beagle_gwas_filtered.vcf.gz
zcat beagle_gwas_filtered.vcf.gz | cut -f1-9,191-200 | gzip > beagle_gwas_filtered.vcf.gz

bcftools index -t beagle_gwas_filtered.vcf.gz

echo
echo "*** Running test analysis with \"gt=\" argument ***"
echo
java –Xmx20g  -jar beagle.27Feb25.75f.jar gt=beagle_gwas_filtered.vcf.gz out=gwas_imputed

bcftools index -t gwas_imputed.vcf.gz
