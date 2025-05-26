#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name scaffold
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --time=48:00:00
#SBATCH -o scaffolding.%j.out
#SBATCH -A lp_svbelleghem

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate thesis 

#YaHs will take the contig sequences (.fa) and the HiC aligned to the contigs (.bam produced in step 4) and scaffold them 
yahs purged.fa mapped.marked.bam  -q 40 -l 10000 -R 2
# q - Increase minimum mapping quality to 40. Removes noisy/ambiguous Hi-C contacts.
# l - Scaffold only contigs > 10 kb. Smaller contigs tend to be repetitive or misassembled and can fragment scaffolds.
# R - Run 2 rounds of scaffolding at each resolution level. May improve robustness without overfitting.

#final scaffolds (.fa) are used to run the ARIMA pipeline again, the mapping will be using the scaffolds as the reference (instead of the contigs)
mv yahs.out_scaffolds_final.fa P_chalceus_REF1.fa

module load SAMtools/1.13-GCC-10.3.0

#finally, the newly re-named assembly is also indexed in order to run the mapping pipeline once more to generate the HiC contact maps
samtools faidx P_chalceus_REF1.fa && cut -f1,2 P_chalceus_REF1.fa.fai > P_chalceus_REF1.fa.genome
