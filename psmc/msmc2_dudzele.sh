#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name 2_msmc2
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=12
#SBATCH --time=72:00:00
#SBATCH --output=2_msmc2.%j.out
#SBATCH -A lp_svbelleghem

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh

module load BCFtools/1.15.1-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load SciPy-bundle/2023.07-gfbf-2023a
conda activate msmc2

#submit with sbatch -a 1-10 (or however many samples in the array)
indID=$((SLURM_ARRAY_TASK_ID -1))

REF=sorted_prim_dud.fasta
IN=/scratch/leuven/357/vsc35707/psmc
OUT=/scratch/leuven/357/vsc35707/psmc/msmc
OUT2=/scratch/leuven/357/vsc35707/psmc/msmc-final

# Sample IDs
SAMPLE=(GC129388 GC129395 GC129400 GC129406 GC129413 GC129417 GC136116 GC136110 GC136117 GC136123)
#belgium SW/LW, france, portugal, spain, heist (2000/2018)

#bam files need to be indexed before starting
samtools index -M $IN/$(echo "${SAMPLE[indID]}").dudPrim.filtered.sorted.nd.bam $IN/$(echo "${SAMPLE[indID]}").dudPrim.filtered.sorted.nd.bam.bai

#Chromosome names
for SCAF in CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag

do
bcftools mpileup -q 20 -Q 20 -C 50 -r $SCAF -f $REF $IN/$(echo "${SAMPLE[indID]}").dudPrim.filtered.sorted.nd.bam | bcftools call -c - | bamCaller.py 30 $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.mask.bed.gz | gzip -c > $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.vcf.gz
generate_multihetsep.py --mask=$OUT/$(echo "${SAMPLE[indID]}")_$SCAF.mask.bed.gz $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.vcf.gz > $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.txt
done

for SCAF in CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag
do
COMMAND="$COMMAND $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.txt"
done

#this script was written for msmc v2.1.4 
msmc2_Linux -t 12 -o $OUT2/$(echo "${SAMPLE[indID]}") $COMMAND

#Sample names
for FINAL in GC129388 GC129395 GC129400 GC129406 GC129413 GC129417 GC136116 GC136110 GC136117 GC136123
do
python MSMC_plotInput.py -I $OUT2/$FINAL.final.txt -u 2.1e-09 -g 1 > $OUT2/$FINAL.final.Rin.txt
done

cp $OUT2/$FINAL.final.Rin.txt $VSC_DATA

#script modified from S. Van Belleghem 2024
