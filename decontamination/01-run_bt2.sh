#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name blobtools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=72:00:00
#SBATCH -o blobtools.%j.out
#SBATCH -A lp_svbelleghem

cd /scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud/

conda activate btk

ASSEMBLY=sorted_prim_dud
MAPPED=mapped_RNA_dud
TAXIDLIST=microbial.ids
WORKING_DIR=/scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud/

#this will generate 5 JSON files (gc.json, identifiers.json, length.json, meta.json, and ncount.json)
blobtools create --fasta $ASSEMBLY.fasta $WORKING_DIR

module load BLAST+/2.13.0-gompi-2022a

export BLASTDB=/scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud/nt

#check if BLAST correctly identifies nt as a BLAST database 
blastdbcmd -db nt -info

blastn -db nt \
       -query $ASSEMBLY.fasta \
       -outfmt "6 qseqid staxids bitscore std" \
       -max_target_seqs 1 \
       -max_hsps 1 \
       -evalue 1e-5 \
       -num_threads 32 \
       -taxidlist $TAXIDLIST \
       -out $ASSEMBLY.ncbi.blastn.out
    
blobtools add \
    --hits $ASSEMBLY.ncbi.blastn.out \
    --taxrule bestsumorder \
    --taxdump /scratch/leuven/357/vsc35707/blobtools/sorted_prim_dud/taxdump \
     $WORKING_DIR

blobtools add \
    --cov $MAPPED.sorted.filtered.bam $WORKING_DIR

blobtools add \
    --busco $ASSEMBLY.busco.bacteria.full_summary.tsv \
    --busco $ASSEMBLY.busco.archaea.full_summary.tsv \
    --busco $ASSEMBLY.busco.fungi.full_summary.tsv \
     $WORKING_DIR

blobtools filter \
     --param length--Min=1000 \
     --param bestsumorder_phylum--Keys=no-hit \
     --fasta $ASSEMBLY.fasta \
     --summary STDOUT \
      $WORKING_DIR
