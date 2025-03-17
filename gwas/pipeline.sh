#SBATCH --cluster=genius 
#SBATCH --job-name GWAS 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o GWAS.%j.out

SOURCE_DIR=
POPULATION=
PHENOTYPE=

#LOAD ALL THE NECESSARY PACKAGES
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
module load tabix/0.2.6-GCCcore-6.4.0
module load Python/3.7.0-foss-2018a
module load Java/11.0.20
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

# Pipeline for running a simple GWAS analysis based on WGS/RADtag data and quantitative phenotypes

# Start with a VCF file with called SNPs, make sure the sample names do not contain any underscores ("_") or any other symbols that might interfere with downstream analyses

echo
echo "*** Check missingness per individual ***"
echo

conda activate vcftools

vcftools --gzvcf $POPULATION.vcf.gz --missing-indv

echo
echo "*** Filter your VCF file ***"
echo

vcftools --gzvcf $POPULATION.vcf.gz --max-missing 0.9 --minQ 30 --maf 0.05 --remove-indels --recode --stdout | bgzip > filtered_$POPULATION.vcf.gz
# At least 90% of individuals must have a called genotype
# Minimum quality score of 30
# Minor allele frequency (MAF) threshold of 0.05
# Exclude indels, keeping only SNPs

# Index VCF file
tabix -p vcf filtered_$POPULATION.vcf.gz

echo
echo "*** Split multiallelic SNPs ***"
echo

# Split the multiallelic SNPs into multiple biallelic ones
bcftools norm -m -any -o split_$POPULATION.vcf.gz -Oz filtered_$POPULATION.vcf.gz

echo
echo "*** Imputing missing genotypes ***"
echo

# Java is needed to run BEAGLE (load the appropriate module)
java -Xmx20480m -jar /data/leuven/357/vsc35707/beagle.27Feb25.75f.jar gt=split_$POPULATION.vcf.gz out=imputed_$POPULATION

# Index the imputed vcf file 
bcftools index -t imputed_$POPULATION.vcf.gz

# Gives a brief summary of stats, imputed genotypes 
bcftools stats imputed_$POPULATION.vcf.gz | grep -E "SN|TSTV"
