#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=barbate_fst
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o barbate_fst.%j.out
#SBATCH --array=1-10 

# Load the programs we will use
module load Python/3.7.0-foss-2018a

REFNAME=dudPrim
chrom=$SLURM_ARRAY_TASK_ID  # Assign chromosome number based on job array ID

echo "================="
echo "Processing chromosome $chrom"

pop1=Bar2
pop2=Bar4

python popgenWindows_egglib.py -w 50000 -s 50000 --minSites 1000 --maxMissing 0.25 \
-T 10 --windType coordinate -f phased \
-g Pogonus_Barbate_$REFNAME.chr_$chrom.H.calls.gz \
--popsFile Pogonus_pops_Barbate.txt \
-o Pogonus_Barbate_$REFNAME.chr_$chrom.stats_${pop1}_${pop2}_w50000_s50000_eggStats.stats \
-p $pop1 \
-p $pop2 \
-eggB FstWC,Dxy -eggW S,Pi,thetaW,D
