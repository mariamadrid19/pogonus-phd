#!/bin/bash -l 
#SBATCH --cluster=genius 
#SBATCH --job-name barbate_fst
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=20 
#SBATCH --time=12:00:00 
#SBATCH -A lp_svbelleghem
#SBATCH -o barbate_fst.%j.out
#SBATCH --array=1-11

# This variable will store the job array number minus 1, so we can use it to get a sample from the samples list (index starts at 0)
ID=$((SLURM_ARRAY_TASK_ID -1))

# Load the programs we will use
module load Python/3.7.0-foss-2018a
#module load tabix
export PYTHONPATH=$PYTHONPATH:/vsc-hard-mounts/leuven-data/350/vsc35085/programs/genomics_general

REFNAME=dudPrim
chrom=1

echo "================="

pop1=Bar2
pop2=Bar4

python popgenWindows_egglib.py -w 50000 -s 50000 --minSites 1000 --maxMissing 0.25 \
-T 10 --windType coordinate -f phased \
-g Pogonus_Barbate_$REFNAME.chr_$chrom.H.calls.gz \
--popsFile Pogonus_pops_Barbate.txt \
-o Pogonus_Barbate_$REFNAME.chr_$chrom.stats_$(echo "${pop1[ID]}")_$(echo "${pop2[ID]}")_w50000_s50000_eggStats.stats \
-p $(echo "${pop1[ID]}") \
-p $(echo "${pop2[ID]}") \
-eggB FstWC,Dxy -eggW S,Pi,thetaW,D
