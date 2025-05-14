#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name purging_Plit
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --time=24:00:00
#SBATCH -o purging_Plit.%j.out
#SBATCH -A lp_svbelleghem

conda deactivate
compleasm.py run -a Pogonus_litt_l3.asm.hic.p_ctg.fa -o results_l3/ -l coleoptera -t 32

conda activate thesis
module load matplotlib/3.7.0-gfbf-2022b

#assembly and reads
PRIMASSEMBLY=Pogonus_litt_l3.asm.hic.p_ctg.fa
FASTQ=bc2042.fastq.gz

#Run minimap2 to align pacbio data and generate paf files, then calculate read depth histogram and base-level read depth
minimap2 -t 36 -xmap-hifi $PRIMASSEMBLY $FASTQ | gzip -c - > $PRIMASSEMBLY.paf.gz
pbcstat $PRIMASSEMBLY.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log
hist_plot.py -c cutoffs PB.stat hist.png

#Split an assembly and do a self-self alignment
split_fa $PRIMASSEMBLY > $PRIMASSEMBLY.split
minimap2 -t 36 -xasm5 -DP $PRIMASSEMBLY.split $PRIMASSEMBLY.split | pigz > $PRIMASSEMBLY.split.self.paf.gz

#Purge haplotigs and overlaps
purge_dups -2 -T cutoffs_min15 -c PB.base.cov $PRIMASSEMBLY.split.self.paf.gz > dups.bed 2> purge_dups.log

#Get purged primary and haplotig sequences from draft assembly
get_seqs -e dups.bed $PRIMASSEMBLY

compleasm.py run -a purged.fa -o results_purged/ -l coleoptera -t 32
