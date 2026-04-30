awk -v BIN=100000 '
BEGIN { OFS="\t" }

$0 !~ /^#/ && $3 == "exon" {
    chrom = $1
    exon_start = $4 - 1   # convert GTF 1-based start to BED-like 0-based
    exon_end = $5

    for (bin_start = int(exon_start / BIN) * BIN; bin_start < exon_end; bin_start += BIN) {
        bin_end = bin_start + BIN

        overlap_start = exon_start > bin_start ? exon_start : bin_start
        overlap_end = exon_end < bin_end ? exon_end : bin_end

        overlap = overlap_end - overlap_start

        if (overlap > 0) {
            key = chrom SUBSEP bin_start SUBSEP bin_end
            exon_bp[key] += overlap
            chrom_seen[chrom] = 1
            if (bin_end > chrom_max[chrom]) chrom_max[chrom] = bin_end
        }
    }
}

END {
    for (chrom in chrom_seen) {
        for (bin_start = 0; bin_start < chrom_max[chrom]; bin_start += BIN) {
            bin_end = bin_start + BIN
            key = chrom SUBSEP bin_start SUBSEP bin_end

            percent = 100 * exon_bp[key] / BIN

            print chrom, bin_start, bin_end, percent
        }
    }
}
' braker.gtf | sort -k1,1V -k2,2n > exon_density_100kb.tsv
