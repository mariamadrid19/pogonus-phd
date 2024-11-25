#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name psmc_dud
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=48:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH --array=1-10

ID=$((SLURM_ARRAY_TASK_ID -1))

# Sample IDs
samples=(\
GC129388 GC129394 GC129401 GC129406 \
GC129412 GC129417 GC129424 GC136078 \
GC136096 GC136084 GC136090 GC136105 \
GC136107 GC136109)

REF=sorted_prim_dud.fasta

module load BCFtools/1.15.1-GCC-11.3.0

bcftools mpileup -C50 -f $REF $(echo "${samples[ID]}").dudPrim.filtered.sorted.nd.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > $(echo "${samples[ID]}")_prim_dud.fq.gz

gunzip $(echo "${samples[ID]}")_prim_dud.fq.gz 
 
fq2psmcfa -q 20 -g 100 -s 10 $(echo "${samples[ID]}")_prim_dud.fq > $(echo "${samples[ID]}")_prim_dud.psmcfa

psmc -N25 -t12 -r5 -p "4+25*2+4+6" -o $(echo "${samples[ID]}")_prim_dud.psmc $(echo "${samples[ID]}")_prim_dud.psmcfa
