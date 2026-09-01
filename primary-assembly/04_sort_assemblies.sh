#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name sort_scaffolds
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH --mem=32G
#SBATCH -o contigous_scaffolds.%j.out
#SBATCH -A lp_svbelleghem

source /data/leuven/357/vsc35707/miniconda3/etc/profile.d/conda.sh
conda activate ragtag

cd /scratch/leuven/357/vsc35707/pogonus-genomes

reference_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/chalceus/Pchalceus_LW_final_CHR.fasta"
olivaceus_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/olivaceus/purged.fa"
littoralis_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/littoralis/P_littoralis_REF1.fa"
gilvipes_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/gilvipes/P_gilvipes_REF1.fa"

#First step is to sort the assemblies based on the reference assembly (only chromosomes)
ragtag.py scaffold "$reference_genome" "$olivaceus_genome" -t 24 -o olivaceus_ref
ragtag.py scaffold "$reference_genome" "$littoralis_genome" -t 24 -o littoralis_ref
ragtag.py scaffold "$reference_genome" "$gilvipes_genome" -t 24 -o gilvipes_ref
