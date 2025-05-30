#! /bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name PretextMap
#SBATCH --nodes=1
#SBATCH --cpus-per-task=36
#SBATCH --time=12:00:00
#SBATCH -o pretext.%j.out
#SBATCH -A lp_svbelleghem

conda activate thesis

# this uses the final bam file produced by the second mapping (to the scaffolded assembly) to produce a PRETEXT map that will be visualized in PretextView
samtools view -h mapped.marked.bam | PretextMap -o P_chalceus_REF1.pretext --sortby length --mapq 10 --highRes

cp P_chalceus_REF1.pretext $VSC_DATA
