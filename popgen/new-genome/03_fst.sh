#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=fst
#SBATCH --ntasks-per-node=20
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o fst.%A_%a.out
#SBATCH --array=0-54

# Array ID breakdown
CHR_INDEX=$((SLURM_ARRAY_TASK_ID / 5 + 1))   # Chromosome 1 to 11
PAIR_INDEX=$((SLURM_ARRAY_TASK_ID % 5))      # Comparison 0 to 4

echo "==================="
echo "Processing CHR${CHR_INDEX}, POPULATION COMPARISON ${PAIR_INDEX}"

# Define population pairs
pop1=("B_NIE_T" "F_GUE_T" "P_AVE_T" "S_HUE_T" "E_SEV_T")
pop2=("B_DUD_S" "F_GUE_S" "P_AVE_S" "S_COT_S" "F_CAM_S")
POPSFILE=popfileALL.txt

# Paths
CALLS_DIR="final-calls"
OUTDIR="final-stats"
mkdir -p "$OUTDIR"

CALLS_FILE="${CALLS_DIR}/Pchal_Bar_SW.chr_${CHR_INDEX}.H.calls.gz"
OUTFILE="${OUTDIR}/Pchal_Bar_SW.chr_${CHR_INDEX}.pop${PAIR_INDEX}.stats"

# Run popgen analysis
/data/leuven/357/vsc35707/miniconda3/bin/python popgenWindows_egglib.py \
    -w 25000 -s 12500 --minSites 250 --maxMissing 0.25 \
    -T 10 --windType coordinate -f phased \
    -g "$CALLS_FILE" \
    --popsFile "$POPSFILE" \
    -o "$OUTFILE" \
    -p "${pop1[$PAIR_INDEX]}" \
    -p "${pop2[$PAIR_INDEX]}" \
    -eggB FstWC,Dxy -eggW Pi,D

echo "Done with CHR${CHR_INDEX}, comparison: ${pop1[$PAIR_INDEX]} vs ${pop2[$PAIR_INDEX]}"
