echo -e "chrom\tpos\tcM_m\tcM_f\tLG" > linkage_table.tsv

for X in {1..11}; do
    awk -v lg="$X" '
        BEGIN { OFS="\t" }
        NR==FNR {
            chrom[FNR-1] = $1
            pos[FNR-1]   = $2
            next
        }
        !/^#/ {
            idx = $1
            if (idx in chrom) {
                printf "%s\t%s\t%.3f\t%.3f\tLG%02d\n", chrom[idx], pos[idx], $2, $3, lg
            }
        }
    ' snps_t05.txt "map12_chr${X}_mrecom0.txt" >> linkage_table.tsv
done
