module load SAMtools/1.16.1-GCC-11.3.0

module load BCFtools/1.15.1-GCC-11.3.0

samtools mpileup -C50 -uf sorted_prim_dud.fasta GC136107.dudPrim.filtered.sorted.nd.bam | bcftools view -c - | vcfutils.pl vcf2fq -d 10 -D 100 | gzip > diploid.fq.gz
