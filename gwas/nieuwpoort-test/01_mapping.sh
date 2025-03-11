#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name call_SNPs 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o call_snps.%j.out
#SBATCH --array=1-11

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index  starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
#module load BCFtools/1.9-foss-2018a
module load Python/3.7.0-foss-2018a
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

echo "================="

chrom=(CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag)
names=(1 10 2 3 4 5 6 7 8 9)

# Sample IDs (24 samples)
samples=(GC129388 GC129389 GC129390 GC129391 GC129392 GC129393 \
Na_034 PcNP_034 PcNP_037 PcNP_040 PcNP_044 PcNP_045 PcNP_033 \
Nb_002 Nb_006 Nb_007 Nb_015 Nb_019 Nb_025 Nb_031 Nb_033 Nb_038 Nb_062 Nb_095)

# Some folder and file paths to use later
REF=/scratch/leuven/357/vsc35707/GWAS/sorted_prim_dud.fasta
REFNAME=primDud

cd /scratch/leuven/357/vsc35707/GWAS/bams

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
bcftools mpileup -Oz --threads 20 -f $REF $(echo $command) -r $(echo "${chrom[ID]}") | \
bcftools call -m -Oz -o /scratch/leuven/357/vsc35707/GWAS/Pogonus_ALL_$REFNAME.chr_$(echo "${names[ID]}").vcf.gz
