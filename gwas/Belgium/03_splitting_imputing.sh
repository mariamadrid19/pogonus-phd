#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name impute 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o impute.%j.out

# tabix is loaded to use bgzip
module load tabix/0.2.6-GCCcore-6.4.0
# python is loaded to use bcftools
module load Python/3.7.0-foss-2018a
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

POPULATION=Belgium
cd /scratch/leuven/357/vsc35707/GWAS/$POPULATION

echo
echo "*** Splitting multiallelic SNPs ***"
echo

# Split the multiallelic SNPs into multiple biallelic ones
bcftools norm -m -any -o split_variants_$POPULATION.vcf.gz -Oz gwas_filtered_$POPULATION.vcf.gz

echo
echo "*** Imputing missing genotypes with BEAGLE 5.5 ***"
echo

# Java is needed to run BEAGLE (load the appropriate module)
module load Java/11.0.20
java -Xmx20480m -jar beagle.27Feb25.75f.jar gt=split_variants_$POPULATION.vcf.gz out=gwas_imputed_$POPULATION

# Index the imputed vcf file 
bcftools index -t gwas_imputed_$POPULATION.vcf.gz

# Gives a brief summary of stats, imputed genotypes 
bcftools stats gwas_imputed_$POPULATION.vcf.gz | grep -E "SN|TSTV"
