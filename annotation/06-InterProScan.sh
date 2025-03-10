#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name interpro
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o interpro.%j.out
#SBATCH -A lp_svbelleghem

#documentation: https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan

#wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.70-102.0/interproscan-5.70-102.0-64-bit.tar.gz
#wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.70-102.0/interproscan-5.70-102.0-64-bit.tar.gz.md5

# Recommended checksum to confirm the download was successful:
#md5sum -c interproscan-5.70-102.0-64-bit.tar.gz.md5

#tar -pxvzf interproscan-5.70-102.0-*-bit.tar.gz

cd /scratch/leuven/357/vsc35707/genome_annotation/braker/Augustus/interproscan

#python3 setup.py -f interproscan.properties

module load Java/11.0.20

#need to remove the * from the aa file
#sed -i 's/\*//g' augustus.hints.aa
interproscan.sh -i /scratch/leuven/357/vsc35707/genome_annotation/braker/Augustus/augustus.hints.aa --outfile augustus_protein_interpro

interproscan.sh -i /scratch/leuven/357/vsc35707/genome_annotation/braker/braker.fa --outfile braker_protein_interpro
