#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name interpro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o interpro.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/genome_annotation/braker/Augustus/my_interproscan
#wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.70-102.0/interproscan-5.70-102.0-64-bit.tar.gz
#wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.70-102.0/interproscan-5.70-102.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
#md5sum -c interproscan-5.70-102.0-64-bit.tar.gz.md5

#tar -pxvzf interproscan-5.70-102.0-*-bit.tar.gz

cd interproscan-5.70-102.0/

python3 setup.py -f interproscan.properties

f=(0 2000 4000 6000 8000 10000 12000 14000 16000 18000 20000 22000 24000 26000 28000 30000 32000)
 
fID=$((PBS_ARRAYID -1))
 
interproscan.sh -i augustus.hints.h.$(echo "${f[fID]}").fa
