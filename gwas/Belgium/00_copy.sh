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
