#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name=Spain_fst
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o Spain_fst.%A_%a.out
#SBATCH --array=0-19

# ========== SETUP ==========
REFNAME=Pchal

# Array ID breakdown
CHR_INDEX=$((SLURM_ARRAY_TASK_ID / 2 + 1))  # Chromosome 1–10
PAIR_INDEX=$((SLURM_ARRAY_TASK_ID % 2))     # 0 or 1 for pair 1 or 2

echo "==================="
echo "Processing CHR${CHR_INDEX}, POPULATION COMPARISON ${PAIR_INDEX}"

# Define population pairs
pop1=("B_NIE_T" "S_HUE_T")
pop2=("B_DUD_S" "S_COT_S")
POPSFILE=popfile.txt

# Paths
CALLS_DIR="final-calls"
OUTDIR="final-stats"
mkdir -p "$OUTDIR"

CALLS_FILE="${CALLS_DIR}/Pogonus_${REFNAME}.CHR${CHR_INDEX}.H.calls.gz"
OUTFILE="${OUTDIR}/Pogonus_${REFNAME}.CHR${CHR_INDEX}.pop${PAIR_INDEX}.stats"

# Run popgen analysis
/data/leuven/357/vsc35707/miniconda3/bin/python popgenWindows_egglib.py \
    -w 50000 -s 50000 --minSites 1000 --maxMissing 0.25 \
    -T 10 --windType coordinate -f phased \
    -g "$CALLS_FILE" \
    --popsFile "$POPSFILE" \
    -o "$OUTFILE" \
    -p "${pop1[$PAIR_INDEX]}" \
    -p "${pop2[$PAIR_INDEX]}" \
    -eggB FstWC,Dxy -eggW Pi,D

echo "Done with CHR${CHR_INDEX}, comparison: ${pop1[$PAIR_INDEX]} vs ${pop2[$PAIR_INDEX]}"
