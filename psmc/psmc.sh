#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name psmc 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/psmc

module load BCFtools/1.15.1-GCC-11.3.0

bcftools mpileup -C50 -f sorted_prim_dud.fasta GC136107.dudPrim.filtered.sorted.nd.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > GC136107_prim_dud.fq.gz

utils/fq2psmcfa -q20 GC136107_prim_dud.fq.gz > GC136107_prim_dud.psmcfa

psmc -N25 -t12 -r5 -p "4+25*2+4+6" -o GC136107_prim_dud.psmc GC136107_prim_dud.psmcfa

utils/psmc_plot.pl GC136107_prim_dud GC136107_prim_dud.psmc
