#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=barbate_fst
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o barbate_fst.%j.out
#SBATCH --array=1-263

REFNAME=P_chalceus_REF1
chrom=$SLURM_ARRAY_TASK_ID  # Assign chromosome number based on job array ID

echo "================="
echo "Processing chromosome $chrom"

pop1=LW
pop2=SW

/data/leuven/357/vsc35707/miniconda3/bin/python popgenWindows_egglib.py -w 50000 -s 50000 --minSites 1000 --maxMissing 0.25 \
-T 10 --windType coordinate -f phased \
-g Pogonus_$REFNAME.chr_$chrom.H.calls.gz \
--popsFile Pogonus_pops.txt \
-o Pogonus_$REFNAME.chr_$chrom.stats_${pop1}_${pop2}.stats \
-p $pop1 \
-p $pop2 \
-eggB FstWC
