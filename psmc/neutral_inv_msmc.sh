#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name msmc_inv
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --time=24:00:00
#SBATCH --output=inv_msmc.%j.out
#SBATCH --array=1-101
#SBATCH -A lp_edu_eeg_2024

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh

module load BCFtools/1.15.1-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load SciPy-bundle/2023.07-gfbf-2023a
conda activate msmc2

indID=$((SLURM_ARRAY_TASK_ID -1))

REF=sorted_prim_dud.fasta
IN=/scratch/leuven/357/vsc35707/psmc/bam-files
OUT=/scratch/leuven/357/vsc35707/psmc/msmc
OUT2=/scratch/leuven/357/vsc35707/psmc/msmc-final
INV_OUT=/scratch/leuven/357/vsc35707/psmc/inverted-msmc

SAMPLE=(GC129388 GC129389 GC129390 GC129391 GC129392 GC129393 GC129394 GC129395 GC129396 GC129397 GC129398 GC129399 GC129400 GC129401 GC129402 GC129403 GC129404 GC129405 GC129406 GC129407 GC129408 GC129409 GC129410 GC129411 GC129412 GC129413 GC129414 GC129415 GC129416 GC129417 \
GC129418 GC129419 GC129420 GC129421 GC129422 GC129423 GC129424 GC129425 GC129426 GC129427 GC129428 GC129429 GC129430 GC129431 GC129432 GC129433 GC129434 GC129435 GC129437 GC129438 GC129439 GC129440 GC136078 GC136079 GC136080 GC136081 GC136082 GC136083 GC136084 GC136085 GC136086 \
GC136087 GC136088 GC136089 GC136090 GC136091 GC136092 GC136093 GC136094 GC136095 GC136096 GC136097 GC136098 GC136099 GC136100 GC136101 GC136102 GC136103 GC136104 GC136105 GC136106 GC136107 GC136108 GC136109 GC136110 GC136111 GC136112 GC136113 GC136114 GC136115 GC136116 GC136117 \
GC136118 GC136119 GC136120 GC136121 GC136122 GC136123 GC136124 GC136125 GC136126 GC136127 GC136128)


#bam files need to be indexed before starting
samtools index -M $IN/$(echo "${SAMPLE[indID]}").dudPrim.filtered.sorted.nd.bam $IN/$(echo "${SAMPLE[indID]}").dudPrim.filtered.sorted.nd.bam.bai

# Define inverted regions
declare -A inverted_regions
inverted_regions["CM008230.1_RagTag"]="3075068-117176353"
inverted_regions["CM008231.1_RagTag"]="23764-14474522"
inverted_regions["CM008233.1_RagTag"]="676248-71874509"
inverted_regions["CM008234.1_RagTag"]="5374868-45974700"
inverted_regions["CM008235.1_RagTag"]="74493-33562646"
inverted_regions["CM008236.1_RagTag"]="1126238-49965318"
inverted_regions["CM008237.1_RagTag"]="11675768-53227168"
inverted_regions["CM008238.1_RagTag"]="8575148-35074000"
inverted_regions["CM008239.1_RagTag"]="7525414-18420609"
inverted_regions["CM008240.1_RagTag"]="9624908-21575397"

# Function to check if a region is inverted
is_inverted() {
    local scaffold=$1
    local start=$2
    local end=$3
    local inv_start=${inverted_regions[$scaffold]%-*}
    local inv_end=${inverted_regions[$scaffold]#*-}
    if (( start >= inv_start && end <= inv_end )); then
        echo "true"
    else
        echo "false"
    fi
}

# Loop through each scaffold (SCAF)
for SCAF in CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag

do

bcftools mpileup -q 20 -Q 20 -C 50 -r $SCAF -f $REF $IN/$(echo "${SAMPLE[indID]}").dudPrim.filtered.sorted.nd.bam | bcftools call -c - | bamCaller.py 30 $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.mask.bed.gz | gzip -c > $OUT/$(echo "${SAMPLE[indID]}")_$SCAF.vcf.gz

# Generate multi-heterozygosity separation for inverted and neutral regions
    if [ $(is_inverted $SCAF $start $end) == "true" ]; then
        generate_multihetsep.py --mask=$OUT/"${SAMPLE[indID]}_$SCAF.mask.bed.gz" $OUT/"${SAMPLE[indID]}_$SCAF.vcf.gz" > $INV_OUT/"${SAMPLE[indID]}_$SCAF.inverted.txt"
    else
        generate_multihetsep.py --mask=$OUT/"${SAMPLE[indID]}_$SCAF.mask.bed.gz" $OUT/"${SAMPLE[indID]}_$SCAF.vcf.gz" > $INV_OUT/"${SAMPLE[indID]}_$SCAF.neutral.txt"
    fi
done

# Collect all generated text file paths for further processing
COMMAND_INVERTED=""
COMMAND_NON_INVERTED=""
for SCAF in CM008230.1_RagTag CM008231.1_RagTag CM008233.1_RagTag CM008234.1_RagTag CM008235.1_RagTag CM008236.1_RagTag CM008237.1_RagTag CM008238.1_RagTag CM008239.1_RagTag CM008240.1_RagTag
do
    if [ $(is_inverted $SCAF $start $end) == "true" ]; then
        COMMAND_INVERTED="$COMMAND_INVERTED $INV_OUT/${SAMPLE[indID]}_$SCAF.inverted.txt"
    else
        COMMAND_NON_INVERTED="$COMMAND_NON_INVERTED $INV_OUT/${SAMPLE[indID]}_$SCAF.neutral.txt"
    fi
done

# Run MSMC2 for inverted regions
msmc2_Linux -t 12 -o $INV_OUT/$(echo "${SAMPLE[indID]}_inverted") $COMMAND_INVERTED

# Run MSMC2 for neutral regions
msmc2_Linux -t 12 -o $INV_OUT/$(echo "${SAMPLE[indID]}_neutral") $COMMAND_NON_INVERTED

# Generate plot input files for inverted regions
for FINAL in "${SAMPLE[@]}"
do
    if [ -f $INV_OUT/${FINAL}_inverted.final.txt ]; then
        python MSMC_plotInput.py -I $INV_OUT/${FINAL}_inverted.final.txt -u 2.1e-09 -g 1 > $INV_OUT/${FINAL}_inverted.final.Rin.txt
    fi
done

# Generate plot input files for neutral regions
for FINAL in "${SAMPLE[@]}"
do
    if [ -f $INV_OUT/${FINAL}_neutral.final.txt ]; then
        python MSMC_plotInput.py -I $INV_OUT/${FINAL}_neutral.final.txt -u 2.1e-09 -g 1 > $INV_OUT/${FINAL}_neutral.final.Rin.txt
    fi
done
