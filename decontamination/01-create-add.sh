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

blastn -db nt \
       -query sorted_prim_dud.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 10 \
       -max_hsps 1 \
       -evalue 1e-25 \
       -num_threads 16 \
       -out sorted_prim_dud.ncbi.blastn.out.txt
       
blobtools add \
    --hits sorted_prim_dud.ncbi.blastn.out.txt \
    --taxrule bestsumorder \
    --taxdump /scratch/leuven/357/vsc35707/blobtools/taxdump

blobtools add --cov Dud_mapped_IsoSeq.filtered.bam --threads 24 --replace

blobtools create -i sorted_prim_dud.fasta -b mapping_sorted.bam -t sorted_prim_dud.ncbi.blastn.out.txt -o filename
