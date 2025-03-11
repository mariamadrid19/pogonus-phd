#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name merge_Np
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o merge_Np.%j.out

module load Python/3.7.0-foss-2018a
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

cd /scratch/leuven/357/vsc35707/GWAS/Nieuwpoort

bcftools concat -Oz -o merged_Npt.vcf.gz *.vcf.gz
