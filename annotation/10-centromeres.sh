#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name quartet
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=24:00:00
#SBATCH -o quartet.%j.out
#SBATCH -A lp_svbelleghem

mamba activate quartet

GENOME_LW="/scratch/leuven/357/vsc35707/centromeres/genomes/Pchalceus_LW_final.fasta"
REPEATS_LW="/scratch/leuven/357/vsc35707/centromeres/repeats/repeats_LW.gff3"
GENES_LW="/scratch/leuven/357/vsc35707/centromeres/genes/braker_LW.gff3"

GENOME_SW="/scratch/leuven/357/vsc35707/centromeres/genomes/Pchalceus_SW_LW_chromosome_names.fasta"
REPEATS_SW="/scratch/leuven/357/vsc35707/centromeres/repeats/repeats_SW.gff3"
GENES_SW="/scratch/leuven/357/vsc35707/centromeres/genes/braker_SW.gff3"

quartet CentroMiner -i "$GENOME_LW" --TE "$REPEATS_LW" --gene "$GENES_LW" -t 24 -p Pchalceus_LW

quartet CentroMiner -i "$GENOME_SW" --TE "$REPEATS_SW" --gene "$GENES_SW" -t 24 -p Pchalceus_SW
