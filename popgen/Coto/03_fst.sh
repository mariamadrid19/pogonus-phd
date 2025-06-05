#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name=Spain_fst
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o Spain_fst.%j.out
#SBATCH --array=1-263

REFNAME=REF1
scaffold=$SLURM_ARRAY_TASK_ID  # Assign scaffold number based on job array ID

echo "================="
echo "Processing scaffold $scaffold"

pop1=LW
pop2=SW
POPSFILE=popfile.txt

/data/leuven/357/vsc35707/miniconda3/bin/python popgenWindows_egglib.py -w 50000 -s 50000 --minSites 1000 --maxMissing 0.25 \
-T 10 --windType coordinate -f phased \
-g Pogonus_$REFNAME.chr_$scaffold.H.calls.gz \
--popsFile $POPSFILE \
-o Pogonus_$REFNAME.chr_$scaffold.stats \
-p $pop1 \
-p $pop2 \
-eggB FstWC,Dxy -eggW Pi,D

echo "Done processing scaffold $scaffold"
