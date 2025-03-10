#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name impute 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o impute.%j.out

# Beagle 5.5
zcat gwas_filtered.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > beagle_gwas_filtered.vcf.gz
zcat beagle_gwas_filtered.vcf.gz | cut -f1-9,191-200 | gzip > beagle_gwas_filtered.vcf.gz

echo
echo "*** Running test analysis with \"gt=\" argument ***"
echo
java –Xmx20g  -jar beagle.27Feb25.75f.jar gt=beagle_gwas_filtered.vcf.gz out=beagle_gwas.gt
