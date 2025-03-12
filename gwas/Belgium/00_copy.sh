#!/bin/bash -l
#SBATCH --cluster=genius 
#SBATCH --job-name copy_BE
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=4:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o copy_BE.%j.out

POPULATION=Belgium
cd /scratch/leuven/357/vsc35707/GWAS/$POPULATION

# Define the directory containing the BAM files
SOURCE_DIR="/scratch/leuven/357/vsc35707/GWAS/bams"

# Loop through each sample name in sample_names.txt
while read sample; do
  # Define file paths
  BAM_FILE="$sample.primDud.filtered.sorted.nd.bam"
  BAI_FILE="$BAM_FILE.bai"

  # Check if the BAM file is already in the current directory, if not, copy it
  if [[ ! -f "$BAM_FILE" ]]; then
    cp "$SOURCE_DIR/$BAM_FILE" . 2>/dev/null && echo "Copied $BAM_FILE"
  else
    echo "$BAM_FILE already exists, skipping."
  fi

  # Check if the BAI file is already in the current directory, if not, copy it
  if [[ ! -f "$BAI_FILE" ]]; then
    cp "$SOURCE_DIR/$BAI_FILE" . 2>/dev/null && echo "Copied $BAI_FILE"
  else
    echo "$BAI_FILE already exists, skipping."
  fi

done < sample_names_BE.txt

echo "File transfer complete."

cd /scratch/leuven/357/vsc35707/GWAS/$POPULATION

bcftools view -s GC129388,GC129389,GC129390,GC129391,GC129392,GC129393,Na_034,PcNP_034,PcNP_037,PcNP_040,PcNP_044,\
PcNP_045,PcNP_033,Nb_002,Nb_006,Nb_007,Nb_015,Nb_019,Nb_025,Nb_031,Nb_033,Nb_038,Nb_062,Nb_095,GC129427,GC129428,\
GC129429,GC129430,GC129431,GC129432,GC129433,GC129434,GC129435,GC129437,GC129438,GC136117,GC136119,GC136120,GC136121,GC136122,GC136123,GC136124,GC136125,GC136126,GC136127,GC136128 \
-Oz -o Pogonus_Np_Hst.vcf.gz /scratch/leuven/357/vsc35707/GWAS/merged_variants.vcf.gz

bcftools view -s GC129388,GC129389,GC129390,GC129391,GC129392,GC129393,Na_034,PcNP_034,PcNP_037,PcNP_040,PcNP_044,PcNP_045,PcNP_033,\
Nb_002,Nb_006,Nb_007,Nb_015,Nb_019,Nb_025,Nb_031,Nb_033,Nb_038,Nb_062,Nb_095,\
GC129394,GC129395,GC129396,GC129397,GC129398,GC129399,GC136084,GC136085,GC136086,GC136087,GC136088,GC136089,\
Db_101,Db_102,Db_103,Pc_DZ_001,Pc_DZ_008,Pc_DZ_012,Pc_DZ_016,Pc_DZ_018,Pc_DZ_022,Pc_DZ_024,Pc_DZ_027,Pc_DZ_030,\
GC129427,GC129428,GC129429,GC129430,GC129431,GC129432,GC129433,GC129434,GC129435,GC129437,GC129438,GC136117,GC136119,GC136120,GC136121,GC136122,GC136123,GC136124,GC136125,GC136126,GC136127,GC136128 \
-O z -o Pogonus_BE.vcf.gz merged_variants.vcf.gz


bcftools query -l Pogonus_Np_Hst.vcf.gz
bcftools index Pogonus_Np_Hst.vcf.gz

bcftools query -l Pogonus_BE.vcf.gz
bcftools index Pogonus_BE.vcf.gz
