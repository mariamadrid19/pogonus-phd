#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name download_db
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o download_db.%j.out
#SBATCH -A lp_svbelleghem

conda activate btk

module load BLAST+/2.13.0-gompi-2022a

blastn -query sorted_prim_dud.fasta -db /scratch/leuven/357/vsc35707/blobtools/nt -outfmt "6 qseqid staxids bitscore std" -max_target_seqs 1 -max_hsps 1 -num_threads 4 -evalue 1e-25 -out sorted_prim_dud.ncbi.blastn.out.txt

blobtools add \
    --hits sorted_prim_dud.ncbi.blastn.out.txt \
    --taxrule bestsumorder \
    --taxdump /scratch/leuven/357/vsc35707/blobtools/taxdump

blobtools create -i sorted_prim_dud.fasta -b mapping_sorted.bam -t sorted_prim_dud.ncbi.blastn.out.txt -o filename
