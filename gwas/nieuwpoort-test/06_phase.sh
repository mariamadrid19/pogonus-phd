#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name phasing 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o phasing.%j.out

module load BCFtools/1.12-GCC-10.3.0

bcftools view -r CHR3 P_chalceus_NP25_BarSW_merged_filtered_mm80_multiSplit.vcf.gz -o chr3_multi_split.vcf.gz -Oz

tabix -p vcf chr3_multi_split.vcf.gz

# Phasing with beagle, no imputing
module load Java/11.0.20
java -Xmx20480m -jar /data/leuven/357/vsc35707/beagle.27Feb25.75f.jar \
  gt=chr3_multi_split.vcf.gz \
  out=chr3_phased \
  nthreads=16 \
  impute=false

tabix -p vcf chr3_phased.vcf.gz

pos_min=38310157
pos_max=38330157

bcftools view -r CHR3:$pos_min-$pos_max chr3_phased.vcf.gz -Oz -o chr3_phased_inversion.vcf.gz

tabix -p vcf chr3_phased_inversion.vcf.gz
