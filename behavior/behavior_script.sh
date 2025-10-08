#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name NAME
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH -A lp_svbelleghem
#SBATCH -o NAME.%j.out

CONVERT FASTQ TO FASTA
cat FILE.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > FILE.fasta

ADD SOMETHING TO PATH
export PATH="/scratch/leuven/357/vsc35707/directory/where/the/program/is:$PATH"
export PATH="/scratch/leuven/357/vsc35707/dudzele_pogonus/gfastats/build/bin:$PATH"

COPY FROM CLUSTER TO COMPUTER
scp -r vsc35707@login-genius.hpc.kuleuven.be:/scratch/leuven/357/vsc35707/ /Users/mariamadrid/Downloads/

COPY FROM COMPUTER TO CLUSTER
scp fileNAME vsc35707@login-genius.hpc.kuleuven.be:/scratch/leuven/357/vsc35707/

SUBMIT INTERACTIVE JOB
srun -M genius -A lp_svbelleghem -n1 -c1 -t60 --mem=128G --pty bash
