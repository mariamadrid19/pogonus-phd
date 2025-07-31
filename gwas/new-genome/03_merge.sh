#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name merge 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o merge_vcfs.%j.out

module load BCFtools/1.12-GCC-10.3.0
export BCFTOOLS_PLUGINS=/data/leuven/357/vsc35707/bcftools/plugins

cd /scratch/leuven/357/vsc35707/GWAS/Nieuwpoort/final-vcfs

bcftools concat -Oz -o Pchal_Bar_SW.filtered.thinned.multiSplit.vcf.gz *.vcf.gz
