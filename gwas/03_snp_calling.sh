#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name snps 
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


# Sample IDs (total of 249 samples)
samples=(\
GC129388 GC129389 GC129390 GC129391 GC129392 GC129393 GC129394 GC129395 GC129396 GC129397 \
GC129398 GC129399 GC129400 GC129401 GC129402 GC129403 GC129404 GC129405 GC129406 GC129407 \
GC129408 GC129409 GC129410 GC129411 GC129412 GC129413 GC129414 GC129415 GC129416 GC129417 \
GC129418 GC129419 GC129420 GC129421 GC129422 GC129423 GC129424 GC129425 GC129426 GC129427 \
GC129428 GC129429 GC129430 GC129431 GC129432 GC129433 GC129434 GC129435 GC129437 GC129438 \
GC129439 GC129440 GC136078 GC136079 GC136080 GC136081 GC136082 GC136083 GC136084 GC136085 \
GC136086 GC136087 GC136088 GC136089 GC136090 GC136091 GC136092 GC136093 GC136094 GC136095 \
GC136096 GC136097 GC136098 GC136099 GC136100 GC136101 GC136102 GC136103 GC136104 GC136105 \
GC136106 GC136107 GC136108 GC136109 GC136110 GC136111 GC136112 GC136113 GC136114 GC136115 \
GC136116 GC136117 GC136118 GC136119 GC136120 GC136121 GC136122 GC136123 GC136124 GC136125 \
GC136126 GC136127 GC136128 Da_022 Da_029 Db_094 Db_098 Db_100 Db_101 Db_102 Db_103 GC3b_071 GC3b_072 GC3b_073 GC3b_074 GC3b_075 GC3b_076 GC3b_077 GC3b_078 GC3b_079 \
GC3b_080 GC3b_081 GC3b_082 GC3b_083 GC3b_084 GC3b_085 GC3b_086 GC3b_087 GC3b_088 GC3b_089 GC3b_090 GC3b_091 GC3b_092 GC3b_093 GC3b_098 GP3_055 \
GP3_059 GP3_067 GP3_068 GP3_078 GP3_088 GP3b_047 GP3b_048 GP3b_055 GP3b_056 GP3b_057 GP3b_062 GP3b_063 GP3b_064 GP3b_066 GP3b_067 GP3b_068 GP3b_069 \
GP3b_076 GP3b_077 GP3b_080 GP3b_082 GP3b_084 GP3b_085 Na_034 Nb_002 Nb_006 Nb_007 Nb_015 Nb_019 Nb_025 Nb_031 Nb_033 Nb_038 Nb_062 Nb_095 Pc2311 \
Pc2314 Pc2316 Pc2317 Pc2318 Pc2321 Pc2323 Pc2324 Pc2325 Pc2331 Pc2332 Pc2336 Pc2337 Pc2338 Pc2339 Pc2340 Pc277 Pc297 Pc298 Pc299 Pc300 Pc316 Pc317 \
Pc318 Pc786 Pc787 Pc788 Pc801 Pc805 Pc821 Pc822 Pc823 PcAVE_003 PcAVE_004 PcAVE_005 PcAVE_006 PcAVE_018 PcAVE_019 PcAVE_025 PcAVE_026 PcAVE_031 PcAVE_032 \
PcAVE_033 PcAVE_034 PcAVE_035 PcAVE_036 PcAVE_037 PcAVE_038 PcNP_033 PcNP_034 PcNP_035 PcNP_036 PcNP_037 PcNP_038 PcNP_040 PcNP_041 PcNP_043 PcNP_044 PcNP_045 \
PcNP_046 Pc_DZ_001 Pc_DZ_002 Pc_DZ_004 Pc_DZ_005 Pc_DZ_007 Pc_DZ_008 Pc_DZ_009 Pc_DZ_011 Pc_DZ_012 Pc_DZ_013 Pc_DZ_016 Pc_DZ_018 Pc_DZ_022 Pc_DZ_024 Pc_DZ_027 Pc_DZ_030)

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
