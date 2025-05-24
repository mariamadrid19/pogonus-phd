#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name scaffold
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00
#SBATCH --output=scaffold_logs/scaffold.%j.out
#SBATCH --error=scaffold_logs/scaffold.%j.err
#SBATCH -A lp_svbelleghem

#this is to tell my slurm job where conda is, so that it can activate my environment
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate thesis 

samtools view ./deduplicated_files/Pogonus_littoralis_r.bam | awk '{print $5}' | sort | uniq -c

awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next;} {seqlen += length($0)} END {print seqlen}' purged.fa | sort -nr

#YaHs will take the contig sequences (.fa) and the HiC aligned to the contigs (.bam produced in step 4) and scaffold them 
yahs purged.fa ./deduplicated_files/Pogonus_littoralis_r.bam -q 30 -l 10000 -r 0

#final scaffolds (.fa) are used to run the ARIMA pipeline again, the mapping will be using the scaffolds as the reference (instead of the contigs)

mv yahs.out_scaffolds_final.fa P_littoralis_REF1.fa
assembly-stats P_littoralis_REF1.fa

module load SAMtools/1.13-GCC-10.3.0

#finally, the newly re-named assembly is also indexed in order to run the ARIMA pipeline once more, this time mapping to this scaffolding assembly, to generate the HiC contact maps
samtools faidx P_littoralis_REF1.fa && cut -f1,2 P_littoralis_REF1.fa.fai > P_littoralis_REF1.fa.genome
