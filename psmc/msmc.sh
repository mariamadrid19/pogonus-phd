#!/bin/bash
#SBATCH --mem-per-cpu=12gb
#SBATCH --time=48:00:00
#SBATCH --job-name=psmc
#SBATCH --error=psmc
#SBATCH --output=psmc
#SBATCH --ntasks=1
#SBATCH --partition=mpi
 
module load BCFtools/1.15.1-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0

indID=$((SLURM_ARRAY_TASK_ID -1))
 
REF=sorted_prim_dud.fasta
IN=/scratch/leuven/357/vsc35707/psmc
OUT=/scratch/leuven/357/vsc35707/psmc/msmc
OUT2=/scratch/leuven/357/vsc35707/psmc/msmc-final

# Sample IDs
SAMPLE= (GC129388 GC129394 GC129401 GC129406)
 
#COMMAND=""

for SCAF in #add scaffold names here
do


bcftools mpileup -q 20 -Q 20 -C 50 -u -r $SCAF -f $REF $IN/$(echo "${SAMPLE[indID]}").RG.MarkDupl.sorted.bam  | bcftools call -cgI - | bamCaller.py 30 $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.mask.bed.gz | gzip -c > $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.vcf.gz
generate_multihetsep.py --mask=$OUT/$(echo "${SAMPLE[indID]}")_$SCAF.mask.bed.gz $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.vcf.gz > $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.txt
 
done
 
 
#for SCAF in Herato0101 Herato0301 Herato0310 Herato0401 Herato0403 Herato0405 Herato0411 Herato0419 Herato0501 Herato0503 Herato0508 Herato0601 Herato0606 Herato0609 Herato0701 Herato0801 Herato0821 Herato0901 Herato1001 Herato1003 Herato1005 Herato1007 Herato1102 Herato1108 Herato1201 Herato1202 Herato1301 Herato1408 Herato1411 Herato1505 Herato1507 Herato1524 Herato1601 Herato1603 Herato1605 Herato1701 Herato1703 Herato1705 Herato1708 Herato1801 Herato1803 Herato1805 Herato1807 Herato1901 Herato1904 Herato1906 Herato1908 Herato1910 Herato2001
#do
 
#COMMAND="$COMMAND $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.txt"
 
#done
 
#~/work/Programs/msmc-master/build/msmc -t 8 -o $OUT2/$(echo "${SAMPLE[indID]}") $COMMAND
 
 
#for FINAL in NCS_2005 NCS_2012 NCS_2020 NCS_2023 NCS_2025 NCS_2556 BC2115 BC2124 STRI_WOM_5779 STRI_WOM_5780 STRI_WOM_5781 NCS_1179 NCS_1979 NCS_2080 NCS_2211 NCS_2217 NCS_2574 NCS_2581 NCS_2609
#do
#python ~/scratch/scripts/MSMC_plotInput.py -I $OUT2/$FINAL.final.txt -u 2e-09 -g 0.25 > OUT2/$FINAL.final.Rin.txt
#done


#script modified from S. Van Belleghem
