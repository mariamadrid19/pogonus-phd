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

bcftools view -S GC129388.primDud.filtered.sorted.nd.bam GC129389.primDud.filtered.sorted.nd.bam GC129390.primDud.filtered.sorted.nd.bam GC129391.primDud.filtered.sorted.nd.bam \
GC129392.primDud.filtered.sorted.nd.bam GC129393.primDud.filtered.sorted.nd.bam Na_034.primDud.filtered.sorted.nd.bam PcNP_034.primDud.filtered.sorted.nd.bam \
PcNP_037.primDud.filtered.sorted.nd.bam PcNP_040.primDud.filtered.sorted.nd.bam PcNP_044.primDud.filtered.sorted.nd.bam PcNP_045.primDud.filtered.sorted.nd.bam \
PcNP_033.primDud.filtered.sorted.nd.bam Nb_002.primDud.filtered.sorted.nd.bam Nb_006.primDud.filtered.sorted.nd.bam Nb_007.primDud.filtered.sorted.nd.bam \
Nb_015.primDud.filtered.sorted.nd.bam Nb_019.primDud.filtered.sorted.nd.bam Nb_025.primDud.filtered.sorted.nd.bam Nb_031.primDud.filtered.sorted.nd.bam \
Nb_033.primDud.filtered.sorted.nd.bam Nb_038.primDud.filtered.sorted.nd.bam Nb_062.primDud.filtered.sorted.nd.bam Nb_095.primDud.filtered.sorted.nd.bam \
GC129394.primDud.filtered.sorted.nd.bam GC129395.primDud.filtered.sorted.nd.bam GC129396.primDud.filtered.sorted.nd.bam GC129397.primDud.filtered.sorted.nd.bam \
GC129398.primDud.filtered.sorted.nd.bam GC129399.primDud.filtered.sorted.nd.bam GC136084.primDud.filtered.sorted.nd.bam GC136085.primDud.filtered.sorted.nd.bam \
GC136086.primDud.filtered.sorted.nd.bam GC136087.primDud.filtered.sorted.nd.bam GC136088.primDud.filtered.sorted.nd.bam GC136089.primDud.filtered.sorted.nd.bam \
Db_101.primDud.filtered.sorted.nd.bam Db_102.primDud.filtered.sorted.nd.bam Db_103.primDud.filtered.sorted.nd.bam Pc_DZ_001.primDud.filtered.sorted.nd.bam \
Pc_DZ_008.primDud.filtered.sorted.nd.bam Pc_DZ_012.primDud.filtered.sorted.nd.bam Pc_DZ_016.primDud.filtered.sorted.nd.bam Pc_DZ_018.primDud.filtered.sorted.nd.bam \
Pc_DZ_022.primDud.filtered.sorted.nd.bam Pc_DZ_024.primDud.filtered.sorted.nd.bam Pc_DZ_027.primDud.filtered.sorted.nd.bam Pc_DZ_030.primDud.filtered.sorted.nd.bam \
GC129427.primDud.filtered.sorted.nd.bam GC129428.primDud.filtered.sorted.nd.bam GC129429.primDud.filtered.sorted.nd.bam GC129430.primDud.filtered.sorted.nd.bam \
GC129431.primDud.filtered.sorted.nd.bam GC129432.primDud.filtered.sorted.nd.bam GC129433.primDud.filtered.sorted.nd.bam GC129434.primDud.filtered.sorted.nd.bam \
GC129435.primDud.filtered.sorted.nd.bam GC129437.primDud.filtered.sorted.nd.bam GC129438.primDud.filtered.sorted.nd.bam GC136117.primDud.filtered.sorted.nd.bam \
GC136119.primDud.filtered.sorted.nd.bam GC136120.primDud.filtered.sorted.nd.bam GC136121.primDud.filtered.sorted.nd.bam GC136122.primDud.filtered.sorted.nd.bam \
GC136123.primDud.filtered.sorted.nd.bam GC136124.primDud.filtered.sorted.nd.bam GC136125.primDud.filtered.sorted.nd.bam GC136126.primDud.filtered.sorted.nd.bam \
GC136127.primDud.filtered.sorted.nd.bam GC136128.primDud.filtered.sorted.nd.bam \
-Oz -o /scratch/leuven/357/vsc35707/GWAS/Belgium/Belgium.vcf.gz /scratch/leuven/357/vsc35707/GWAS/merged_variants.vcf.gz

bcftools query -l Belgium.vcf.gz
bcftools index Belgium.vcf.gz
