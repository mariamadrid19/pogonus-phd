#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name barbate_snps
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o barbate_snps.%j.out
#SBATCH --array=1-10

cd /scratch/leuven/357/vsc35707/BAR_mapping

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load SAMtools/1.9-GCC-6.4.0-2.28
module load Python/3.7.0-foss-2018a
module load tabix/0.2.6-GCCcore-6.4.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

echo "================="

# Ten chromosomes
chrom=(CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag)
names=(1 10 2 3 4 5 6 7 8 9)

# Sample IDs (20 samples)
samples=(Bar2_01 Bar2_02 Bar2_03 Bar2_04 Bar2_05 Bar2_06 Bar2_07 Bar2_08 Bar2_09 Bar2_10 Bar4_01 Bar4_02 Bar4_03 Bar4_04 Bar4_05 Bar4_06 Bar4_07 Bar4_08 Bar4_09 Bar4_10)

REFNAME=dudPrim

echo "${samples[ID]}"

# make a single list of all the samples that can be used in the samtools command
ALL_LIST=""
for FILE in ${samples[*]}
do
ALL_LIST="$ALL_LIST $FILE".$REFNAME.filtered.sorted.nd.bam""
done
eval command=\$$(echo ALL_LIST)

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
bcftools mpileup -Oz --threads 20 -f $REF /scratch/leuven/357/vsc35707/BAR_mapping/bams/$(echo $command) -r $(echo "${chrom[ID]}") | \
bcftools call -m -Oz -o /scratch/leuven/357/vsc35707/BAR_mapping/Pogonus_Barbate_$REFNAME.chr_$(echo "${names[ID]}").vcf.gz

vcftools --gzvcf Pogonus_Barbate_$REFNAME.chr_$(echo "${names[ID]}").vcf.gz --recode --remove-indels --minQ 30 --max-missing 0.25 --stdout | bgzip > Pogonus_Barbate_$REFNAME.chr_$(echo "${names[ID]}").filtered.vcf.gz

bgzip Pogonus_Barbate_$REFNAME.chr_$(echo "${names[ID]}").filtered.vcf.gz

python parseVCF.py \
--minQual 30 \
--skipIndels \
-i Pogonus_Barbate_$REFNAME.chr_$(echo "${names[ID]}").filtered.vcf.gz  | gzip > Pogonus_Barbate_$REFNAME.chr_$(echo "${names[ID]}").calls.gz

zcat Pogonus_Barbate_$REFNAME.chr_$(echo "${names[ID]}").calls.gz | sed 's/.dudPrim.filtered.sorted.nd.bam//g' | bgzip > Pogonus_Barbate_$REFNAME.chr_$(echo "${names[ID]}").H.calls.gz
