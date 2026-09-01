#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name asynt_map
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=48:00:00
#SBATCH -o asynt_map.%j.out
#SBATCH -A lp_edu_evol_quant_genetics

conda activate thesis

cd /scratch/leuven/357/vsc35707/pogonus-genomes

SW_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/chalceus/Pchalceus_SW_CHR.fasta"
LW_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/chalceus/Pchalceus_LW_final_CHR.fasta"
olivaceus_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/olivaceus_ref/P_olivaceus_CHR.fa"
littoralis_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/littoralis_ref/P_littoralis_CHR.fasta"
gilvipes_genome="/scratch/leuven/357/vsc35707/pogonus-genomes/gilvipes_ref/P_gilvipes_CHR.fasta"

#map assemblies to each other 
minimap2 -t 24 -x asm5 "$SW_genome" "$LW_genome" | gzip > g1_g2.paf.gz
minimap2 -t 24 -x asm5 "$LW_genome" "$olivaceus_genome" | gzip > g2_g3.paf.gz
minimap2 -t 24 -x asm5 "$olivaceus_genome" "$littoralis_genome" | gzip > g3_g4.paf.gz
minimap2 -t 24 -x asm5 "$littoralis_genome" "$gilvipes_genome" | gzip > g4_g5.paf.gz
