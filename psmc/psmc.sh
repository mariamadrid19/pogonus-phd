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

gunzip GC136107_prim_dud.fq.gz 
 
fq2psmcfa -q 20 -g 100 -s 10 GC136107_prim_dud.fq > GC136107_prim_dud.psmcfa

psmc -N25 -t12 -r5 -p "4+25*2+4+6" -o GC136107_prim_dud.psmc GC136107_prim_dud.psmcfa


#after installing gnuplot in a local machine (e.g. brew install gnuplot on a Mac), download the perl script from the Github and the psmc files
#run the following command, on a local machine!
#this will produce an .eps file (can be opened with Illustrator or other image editing software)
perl psmc_plot.pl GC136107_prim_dud GC136107_prim_dud.psmc
