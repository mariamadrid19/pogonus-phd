#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name snps 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=36
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o call_snps.%j.out
#SBATCH --array=1-192

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index  starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

conda activate variant_tools

echo "================="

chrom=(CHR1 CHR2 CHR3 CHR4 CHR5 CHR6 CHR7 CHR8 CHR9 CHR10 CHR11)
names=(01 02 03 04 05 06 07 08 09 10 11)

# Sample IDs (all samples, 192)
samples=()
for i in $(seq -w 1 192); do
  samples+=("Pc25Np${i}")
done

REF=/scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/Pchalceus_SW.sorted.fasta
REFNAME=Pchal_Bar_SW

cd /scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/bams

# make a single list of all the samples that can be used in the samtools command
ALL_LIST=""
for FILE in "${samples[@]}"; do
    ALL_LIST+=" ${FILE}.${REFNAME}.filtered.sorted.bam"
done
 
command="$ALL_LIST"

echo "Reference: $REF"
echo "Chromosome: ${chrom[ID]}"
echo "Command: $command"

# Ensure chromosome ID is not empty
if [[ -z "${chrom[ID]}" ]]; then
    echo "Error: chrom[ID] is empty!"
    exit 1
fi

# Ensure BAM files exist
if [[ -z "$command" ]]; then
    echo "Error: No BAM files specified!"
    exit 1
fi

# run mpileup
bcftools mpileup -Oz --threads 36 -f $REF $command -r ${chrom[ID]} | \
bcftools call -m -Oz -o /scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/vcfs/P_chalceus_NP25_$REFNAME.chr_${names[ID]}.vcf.gz
