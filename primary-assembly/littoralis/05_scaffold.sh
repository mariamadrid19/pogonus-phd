#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name scaffold
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=36
#SBATCH --time=48:00:00
#SBATCH -o scaffold.%j.out
#SBATCH -A lp_svbelleghem

#this is to tell my slurm job where conda is, so that it can activate my environment
source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate thesis 

#YaHs will take the contig sequences (.fa) and the HiC aligned to the contigs (.bam produced in step 4) and scaffold them 
yahs /scratch/leuven/357/vsc35707/littoralis/final_assembly/purged.fa /scratch/leuven/357/vsc35707/littoralis/final_assembly/deduplicated_files/Pogonus_littoralis_r.bam -q 100000 -l 30 -r 50000

#final scaffolds (.fa) are used to run the ARIMA pipeline again, the mapping will be using the scaffolds as the reference (instead of the contigs)

mv yahs.out_scaffolds_final.fa P_littoralis_REF1.fa

module load SAMtools/1.13-GCC-10.3.0

#finally, the newly re-named assembly is also indexed in order to run the ARIMA pipeline once more, this time mapping to this scaffolding assembly
samtools faidx P_littoralis_REF1.fa && cut -f1,2 P_littoralis_REF1.fa.fai > P_littoralis_REF1.fa.genome
