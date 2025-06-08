#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name lepmap
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=16G
#SBATCH --time=72:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o lepmap.%j.out

# Change to the working directory
cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap

# Define directories
LEPMAP="/data/leuven/357/vsc35707/LepMap3"

# Run ParentCall2
zcat post.gz|java -cp $LEPMAP/bin ParentCall2 data=pedigree.txt posteriorFile=- removeNonInformative=1|gzip >data.call.gz
zcat data.call.gz | awk 'NR>=7' | cut -f 1,2 > snps.txt

# Run Filtering2
zcat data.call.gz|java -cp $LEPMAP/bin Filtering2 data=- dataTolerance=0.01|gzip >data_f_t01.call.gz

# Run SeparateChromosomes2 
zcat data_f_t01.call.gz | java -cp $LEPMAP/bin SeparateChromosomes2 data=- lodLimit=5 >map5.txt

# CleanMap (for LepAnchor)
paste snps.txt map5.txt|awk '(NR>1)' >cleanMap.input

# Run OrderMakers2 (per chromosome, as assigned by 
java -cp /data/leuven/357/vsc35707/LepMap3/bin/ OrderMarkers2 map=map5.txt data=data_f.call chromosome=1 recombination1=0 > order1.txt
