#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name msmc_neutral
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH --output=neutral_msmc.%j.out
#SBATCH --array=1-101
#SBATCH -A lp_edu_eeg_2024

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh

module load BCFtools/1.15.1-GCC-11.3.0
module load SAMtools/1.16.1-GCC-11.3.0
module load SciPy-bundle/2023.07-gfbf-2023a
conda activate msmc2

indID=$((SLURM_ARRAY_TASK_ID - 1))

REF=sorted_prim_dud.fasta
IN=/scratch/leuven/357/vsc35707/psmc/bam-files
NEUTRAL_OUT=/scratch/leuven/357/vsc35707/psmc/neutral-msmc
NEUTRAL_OUT2=/scratch/leuven/357/vsc35707/psmc/neutral-msmc-final

# Samples
SAMPLE=(GC129388 GC129389 GC129390 GC129391 GC129392 GC129393 GC129394 GC129395 GC129396 GC129397 GC129398 GC129399 GC129400 GC129401 GC129402 GC129403 GC129404 GC129405 GC129406 GC129407 GC129408 GC129409 GC129410 GC129411 GC129412 GC129413 GC129414 GC129415 GC129416 GC129417 \
GC129418 GC129419 GC129420 GC129421 GC129422 GC129423 GC129424 GC129425 GC129426 GC129427 GC129428 GC129429 GC129430 GC129431 GC129432 GC129433 GC129434 GC129435 GC129437 GC129438 GC129439 GC129440 GC136078 GC136079 GC136080 GC136081 GC136082 GC136083 GC136084 GC136085 GC136086 \
GC136087 GC136088 GC136089 GC136090 GC136091 GC136092 GC136093 GC136094 GC136095 GC136096 GC136097 GC136098 GC136099 GC136100 GC136101 GC136102 GC136103 GC136104 GC136105 GC136106 GC136107 GC136108 GC136109 GC136110 GC136111 GC136112 GC136113 GC136114 GC136115 GC136116 GC136117 \
GC136118 GC136119 GC136120 GC136121 GC136122 GC136123 GC136124 GC136125 GC136126 GC136127 GC136128)

# Ensure output directories exist
mkdir -p $NEUTRAL_OUT $NEUTRAL_OUT2

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

# Define chromosome lengths
declare -A chromosome_lengths
chromosome_lengths["CM008230.1_RagTag"]=122030211 #chromosome 1
chromosome_lengths["CM008231.1_RagTag"]=23453803 #chromosome 10
chromosome_lengths["CM008233.1_RagTag"]=82969690 #chromosome 2
chromosome_lengths["CM008234.1_RagTag"]=46345649 #chromosome 3 
chromosome_lengths["CM008235.1_RagTag"]=33576320 #chromosome 4 
chromosome_lengths["CM008236.1_RagTag"]=49979521 #chromosome 5 
chromosome_lengths["CM008237.1_RagTag"]=54902676 #chromosome 6 
chromosome_lengths["CM008238.1_RagTag"]=37253630 #chromosome 7 
chromosome_lengths["CM008239.1_RagTag"]=18441779 #chromosome 8 
chromosome_lengths["CM008240.1_RagTag"]=22778183 #chromosome 9 

# Initialize COMMAND
COMMAND=""

# Process neutral regions for each chromosome
for CHROM in "${!inverted_regions[@]}"; do
    # Parse inversion region
    INV_START=$(echo "${inverted_regions[$CHROM]}" | cut -d'-' -f1)
    INV_END=$(echo "${inverted_regions[$CHROM]}" | cut -d'-' -f2)
    CHROM_END=${chromosome_lengths[$CHROM]}

    # Neutral region 1: From 1 to (INV_START - 1)
    if (( INV_START > 1 )); then
        REGION="1-$(($INV_START - 1))"
        bcftools mpileup -q 20 -Q 20 -C 50 -r "${CHROM}:${REGION}" -f $REF $IN/${SAMPLE[indID]}.dudPrim.filtered.sorted.nd.bam | \
        bcftools call -c - | \
        bamCaller.py 30 $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.mask.bed.gz | \
        gzip -c > $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.vcf.gz

        generate_multihetsep.py --mask=$NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.mask.bed.gz \
        $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.vcf.gz > $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.txt

        COMMAND="$COMMAND $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.txt"
    fi

    # Neutral region 2: From (INV_END + 1) to CHROM_END
    if (( INV_END < CHROM_END )); then
        REGION="$(($INV_END + 1))-$CHROM_END"
        bcftools mpileup -q 20 -Q 20 -C 50 -r "${CHROM}:${REGION}" -f $REF $IN/${SAMPLE[indID]}.dudPrim.filtered.sorted.nd.bam | \
        bcftools call -c - | \
        bamCaller.py 30 $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.mask.bed.gz | \
        gzip -c > $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.vcf.gz

        generate_multihetsep.py --mask=$NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.mask.bed.gz \
        $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.vcf.gz > $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.txt

        COMMAND="$COMMAND $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_${REGION}.txt"
    fi
done

# Concatenate neutral region files for each chromosome
for CHROM in "${!inverted_regions[@]}"; do
    # Define the output concatenated file for the chromosome
    CONCAT_FILE=$NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_neutral_combined.txt

    # Find all neutral region files for this chromosome and concatenate them
    cat $NEUTRAL_OUT/${SAMPLE[indID]}_${CHROM}_*.txt > $CONCAT_FILE

    # Optionally, add the concatenated file to the COMMAND variable
    COMMAND="$COMMAND $CONCAT_FILE"
done

# Run MSMC2 with the concatenated neutral regions
msmc2_Linux -t 24 -o $NEUTRAL_OUT2/${SAMPLE[indID]} $COMMAND || { echo "Error in MSMC2"; exit 1; }

# Generate MSMC plot inputs for each sample
for FINAL in "${SAMPLE[@]}"; do
    python MSMC_plotInput.py -I $NEUTRAL_OUT2/$FINAL.final.txt -u 2.1e-09 -g 1 > $NEUTRAL_OUT2/$FINAL.final.Rin.txt || { echo "Error in MSMC_plotInput for $FINAL"; exit 1; }
done
