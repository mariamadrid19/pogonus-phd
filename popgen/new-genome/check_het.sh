#!/bin/bash
#SBATCH --cluster=genius
#SBATCH --job-name=check_hets
#SBATCH --output=check_hets_%A_%a.out
#SBATCH --array=0-9
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH -A lp_svbelleghem

# Define list of files
FILES=(
"Pchal_Bar_SW.chr_1.H.calls.gz"
"Pchal_Bar_SW.chr_2.H.calls.gz"
"Pchal_Bar_SW.chr_3.H.calls.gz"
"Pchal_Bar_SW.chr_4.H.calls.gz"
"Pchal_Bar_SW.chr_5.H.calls.gz"
"Pchal_Bar_SW.chr_6.H.calls.gz"
"Pchal_Bar_SW.chr_7.H.calls.gz"
"Pchal_Bar_SW.chr_8.H.calls.gz"
"Pchal_Bar_SW.chr_9.H.calls.gz"
"Pchal_Bar_SW.chr_10.H.calls.gz"
)

FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

echo "Analyzing file: $FILE"

# Count heterozygous genotypes (e.g., A/T, C/G, etc. but NOT A/A or N/N)
zcat "$FILE" | awk '
BEGIN { FS = "\t"; het_count = 0; site_count = 0 }
NR == 1 { next }  # skip header
{
    site_count++
    for (i = 3; i <= NF; i++) {
        split($i, gt, "/")
        if (gt[1] != gt[2] && gt[1] != "N" && gt[2] != "N") {
            het_count++
            next
        }
    }
}
END {
    print "File:", FILENAME
    print "Sites analyzed:", site_count
    print "Sites with at least one heterozygous call:", het_count
}
