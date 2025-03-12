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

# Sample IDs (70 samples)
samples=(GC129388 GC129389 GC129390 GC129391 GC129392 GC129393 \
Na_034 PcNP_034 PcNP_037 PcNP_040 PcNP_044 PcNP_045 PcNP_033 \
Nb_002 Nb_006 Nb_007 Nb_015 Nb_019 Nb_025 Nb_031 Nb_033 Nb_038 Nb_062 Nb_095 \
GC129394 GC129395 GC129396 GC129397 GC129398 GC129399 GC136084 GC136085 GC136086 GC136087 GC136088 GC136089 \
Db_101 Db_102 Db_103 \
Pc_DZ_001 Pc_DZ_008 Pc_DZ_012 Pc_DZ_016 Pc_DZ_018 Pc_DZ_022 Pc_DZ_024 Pc_DZ_027 Pc_DZ_030 \
GC129427 GC129428 GC129429 GC129430 GC129431 GC129432 GC129433 GC129434 GC129435 GC129437 GC129438 \
GC136117 GC136119 GC136120 GC136121 GC136122 GC136123 GC136124 GC136125 GC136126 GC136127 GC136128)

# Some folder and file paths to use later
REF=/scratch/leuven/357/vsc35707/GWAS/sorted_prim_dud.fasta
REFNAME=primDud
POPULATION=Belgium

cd /scratch/leuven/357/vsc35707/GWAS/$POPULATION
# All the bam files should be here

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
bcftools call -m -Oz -o /scratch/leuven/357/vsc35707/GWAS/Pogonus_$POPULATION_$REFNAME.chr_$(echo "${names[ID]}").vcf.gz
