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

# Load Java (needed to run LEPMAP and LEPANCHOR)
module load Java/21.0.2

# Run ParentCall2
zcat post.gz|java -cp $LEPMAP/bin ParentCall2 data=pedigree.txt posteriorFile=- removeNonInformative=1 | gzip > data.call.gz

# Run Filtering2
zcat data.call.gz | java -cp $LEPMAP/bin Filtering2 data=- dataTolerance=0.01 | gzip > data_f_t01.call.gz
zcat data_f_t01.call.gz | awk 'NR>=7' | cut -f 1,2 > snps.txt

# Run SeparateChromosomes2 
zcat data_f_t01.call.gz | java -cp $LEPMAP/bin SeparateChromosomes2 data=- lodLimit=6 > map6.txt

# Run JoinSingles2All
zcat data_f_t01.call.gz | java -cp $LEPMAP/bin JoinSingles2All map=map6.txt data=- lodLimit=5 iterate=2 >map6js.txt

gzip -d -c data_f_t01.call.gz > data_f_t01.call

# Create the input for CleanMap (for LepAnchor)
paste snps.txt map6js.txt|awk '(NR>1)' > cleanMap.input
