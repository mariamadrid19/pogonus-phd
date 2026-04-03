#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=lepmap
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=24:00:00
#SBATCH -o logs/lepmap.%j.out
#SBATCH -A lp_svbelleghem

module load Java/21.0.8
conda activate variant_tools

# Change to the working directory
cd /scratch/leuven/357/vsc35707/linkage-mapping/lepmap

# Define directories
LEPMAP="/data/leuven/357/vsc35707/LepMap3"

# # Run ParentCall2
zcat post.gz | java -cp $LEPMAP/bin ParentCall2 data=pedigree.txt XLimit=2 posteriorFile=- removeNonInformative=1 | gzip > data.call.gz

# # Run Filtering2
zcat data.call.gz | java -cp $LEPMAP/bin Filtering2 data=- dataTolerance=0.01 | gzip > data_f_t01.call.gz

# # Get the snp names to a file
zcat data_f_t01.call.gz | awk 'NR>=7' | cut -f 1,2 > snps.txt

# Try different lod limits
for lod in 4 5 6 7 8 9 10 11 12 15 20 25 30; do
    zcat data_f_t01.call.gz | java -cp $LEPMAP/bin SeparateChromosomes2 data=- lodLimit=$lod > map${lod}.txt
done

# Try different lod limits for JoinSingles
for lod in 4 5 6 7 8; do
    zcat data_f_t01.call.gz | java -cp $LEPMAP/bin JoinSingles2All \
        map=map12.txt \
        data=- \
        lodLimit=$lod \
        iterate=2 \
        > map12_js${lod}.txt
done

#Order markers LOD12, set male recombination as 0 (males are achiasmatic) 
for chr in {1..11}; do
    echo "Running chromosome $chr..."
    zcat data_f_t01.call.gz | java -cp $LEPMAP/bin OrderMarkers2 \
        map=map12.txt \
        data=- \
        recombination1=0 \
        chromosome=$chr \
        > map12_chr${chr}_mrecom0.txt
done
