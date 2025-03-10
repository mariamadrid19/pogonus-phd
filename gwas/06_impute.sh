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

# Cut only the necessary columns of the vcf file 
zcat gwas_filtered.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > beagle_gwas_filtered.vcf.gz
zcat beagle_gwas_filtered.vcf.gz | cut -f1-9,191-200 | gzip > beagle_gwas_filtered.vcf.gz

bcftools index -t beagle_gwas_filtered.vcf.gz

echo
echo "*** Imputing missing genotypes with BEAGLE 5.5 ***"
echo

# Java is needed to run BEAGLE (load the appropriate module)
# Command to run beagle as "beagle" was added to the .bashrc file
# alias beagle='java -Xmx20g -jar /data/leuven/357/vsc35707/beagle.27Feb25.75f.jar'
module load Java/11.0.20
beagle gt=beagle_gwas_filtered.vcf.gz out=gwas_imputed

# Index the imputed vcf file 
bcftools index -t gwas_imputed.vcf.gz
