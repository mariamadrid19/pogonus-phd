#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name psmc_dud
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00 
#SBATCH -A lp_svbelleghem

ID=$((SLURM_ARRAY_TASK_ID -1))

# Sample IDs
samples=(\
GC129404 GC129388 GC129389 GC129390 GC129391 \
GC129392 GC129393 GC129394 GC129395 GC129396 \
GC129397 GC129398 GC129399 GC129400 GC129401 \
GC129402 GC129403 GC129405 GC129406 GC129407 \
GC129408 GC129409 GC129410 GC129411 GC129412 \
GC129413 GC129414 GC129415 GC129416 GC129417 \
GC129418 GC129419 GC129420 GC129421 GC129422 \
GC129423 GC129424 GC129425 GC129426 GC129427 \
GC129428 GC129429 GC129430 GC129431 GC129432 \
GC129434 GC129435 GC129437 GC129438 GC129439 \
GC129440 GC129433)

REF=sorted_prim_dud.fasta

module load BCFtools/1.15.1-GCC-11.3.0

bcftools mpileup -C50 -f $REF $(echo "${samples[ID]}").dudPrim.filtered.sorted.nd.bam | bcftools call -c | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > $(echo "${samples[ID]}")_prim_dud.fq.gz

gunzip $(echo "${samples[ID]}")_prim_dud.fq.gz 
 
fq2psmcfa -q 20 -g 100 -s 10 $(echo "${samples[ID]}")_prim_dud.fq > $(echo "${samples[ID]}")_prim_dud.psmcfa

psmc -N25 -t12 -r5 -p "4+25*2+4+6" -o $(echo "${samples[ID]}")_prim_dud.psmc $(echo "${samples[ID]}")_prim_dud.psmcfa
