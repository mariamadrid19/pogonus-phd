#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name impute 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o impute.%j.out

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
module load tabix/0.2.6-GCCcore-6.4.0
# tabix is loaded to use bgzip
module load Python/3.7.0-foss-2018a
# python is loaded to use bcftools
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

echo
echo "*** Splitting multiallelic SNPs ***"
echo

# Split the multiallelic SNPs into multiple biallelic ones
bcftools norm -m -any -o split_variants_Npt.vcf.gz -Oz gwas_filtered_Npt.vcf.gz

echo
echo "*** Imputing missing genotypes with BEAGLE 5.5 ***"
echo

# Java is needed to run BEAGLE (load the appropriate module)
module load Java/11.0.20
java -Xmx20480m -jar /data/leuven/357/vsc35707/beagle.27Feb25.75f.jar gt=split_variants_Npt.vcf.gz out=gwas_imputed_Npt

# Index the imputed vcf file 
bcftools index -t gwas_imputed_Npt.vcf.gz

# Gives a brief summary of stats, imputed genotypes 
bcftools stats gwas_imputed_Npt.vcf.gz | grep -E "SN|TSTV"
