#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=Spain_fst
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o Spain_fst.%j.out
#SBATCH --array=1-263

# Set reference name and get scaffold ID from array
REFNAME=P_chalceus_REF1
scaffold=$SLURM_ARRAY_TASK_ID  # Scaffold number based on job array index

echo "================="
echo "Processing scaffold $scaffold"

# Population labels and popfile path
pop1=LW
pop2=SW
POPSFILE=popfile.txt

# Input and output directories
CALLS_DIR="final-calls"
OUTDIR="final-stats"
mkdir -p "$OUTDIR"

# Input file path
CALLS_FILE="${CALLS_DIR}/Pogonus_${REFNAME}.scaffold_${scaffold}.H.calls.gz"

# Output prefix
OUTFILE="${OUTDIR}/Pogonus_${REFNAME}.scaffold_${scaffold}.stats"

# Run popgen analysis
/data/leuven/357/vsc35707/miniconda3/bin/python popgenWindows_egglib.py \
    -w 50000 -s 50000 --minSites 1000 --maxMissing 0.25 \
    -T 10 --windType coordinate -f phased \
    -g "$CALLS_FILE" \
    --popsFile "$POPSFILE" \
    -o "$OUTFILE" \
    -p "$pop1" \
    -p "$pop2" \
    -eggB FstWC,Dxy -eggW Pi,D

echo "Done processing scaffold $scaffold"
