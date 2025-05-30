#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name synteny
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=36
#SBATCH --time=72:00:00
#SBATCH -o synteny.%j.out
#SBATCH -A lp_svbelleghem

awk '/^>/ { if (count++ >= 20) exit } { print }' P_chalceus_REF1_sorted.fa  > 20_P_chalceus_REF1.fa 

#We can then study the synteny between them 
conda activate ntsynt

ntSynt 11_old_ref.fa 20_P_chalceus_REF1.fa  -p P_chal_SW -t 36 -d 30

python denovo_synteny_block_stats.py --tsv P_chal_SW.synteny_blocks.tsv --fai 11_old_ref.fa.fai 20_P_chalceus_REF1.fa.fai

python sort_ntsynt_blocks.py --synteny_blocks P_chal_SW.synteny_blocks.tsv --sort_order 11_old_ref.fa.fai 20_P_chalceus_REF1.fa.fai --fais > P_chal_SW.synteny_blocks.sorted.tsv

python format_blocks_gggenomes.py --fai 11_old_ref.fa.fai 20_P_chalceus_REF1.fa.fai --prefix P_chal_SW --blocks P_chal_SW.synteny_blocks.sorted.tsv --length 100 --colour 11_old_ref.fa

cp P_chal_SW.links.tsv $VSC_DATA
cp P_chal_SW.sequence_lengths.tsv $VSC_DATA

# to run it locally 
# Rscript plot_synteny_blocks_gggenomes.R -s P_chal_SW.sequence_lengths.tsv -l P_chal_SW.links.tsv --scale 25000000 --p P_chal_SW
