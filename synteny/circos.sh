#!/bin/bash -l
#SBATCH --clusters=wice
#SBATCH --job-name=circos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=36
#SBATCH --mem=120G
#SBATCH --time=60:00:00
#SBATCH --account=lp_edu_eeg_2026
#SBATCH --output=circos.%j.out

###############################################################################
# Two-genome Circos synteny plot
#
# PAF orientation:
#   Query  = Pchalceus_LW_final.fasta
#   Target = Pchalceus_SW_final.fasta
#
# Link colors:
#   Blue = forward-orientation alignment
#   Red  = reverse-orientation alignment
###############################################################################

# Input files
LW_FASTA="Pchalceus_LW_final.fasta"
SW_FASTA="Pchalceus_SW_final.fasta"
LW_FAI="${LW_FASTA}.fai"
SW_FAI="${SW_FASTA}.fai"
PAF="Pchalceus.paf.gz"

# Output files
KARYOTYPE="karyotype.txt"
LINKS="links.txt"
CONFIG="circos.conf"
OUTPUT_PREFIX="Pchalceus_LW_vs_SW"

# Filtering parameters
MIN_SPAN=50000
MAX_DIVERGENCE=0.05
MIN_MAPQ=20

# Prefixes used to distinguish chromosomes with identical names
QUERY_PREFIX="LW_"
TARGET_PREFIX="SW_"

###############################################################################
# 1. Check required programs and input files
###############################################################################

for program in awk gzip circos; do
    if ! command -v "$program" >/dev/null 2>&1; then
        echo "ERROR: Required program not found: $program" >&2
        exit 1
    fi
done

for input_file in "$LW_FASTA" "$SW_FASTA" "$LW_FAI" "$SW_FAI" "$PAF"; do
    if [[ ! -s "$input_file" ]]; then
        echo "ERROR: Missing or empty input file: $input_file" >&2
        exit 1
    fi
done

echo "All required files and programs were found."

###############################################################################
# 2. Verify the PAF orientation
###############################################################################

# Read the first PAF record.
first_record=$(gzip -cd "$PAF" | sed -n '1p')

if [[ -z "$first_record" ]]; then
    echo "ERROR: The PAF file contains no alignments." >&2
    exit 1
fi

read -r query_name query_length strand target_name target_length < <(
    awk 'BEGIN{FS="\t"}
         {
             print $1, $2, $5, $6, $7
         }' <<< "$first_record"
)

# Obtain the corresponding sequence lengths from the FASTA indexes.
lw_index_length=$(
    awk -v sequence="$query_name" '
        $1 == sequence {
            print $2
            exit
        }' "$LW_FAI"
)

sw_index_length=$(
    awk -v sequence="$target_name" '
        $1 == sequence {
            print $2
            exit
        }' "$SW_FAI"
)

echo
echo "First PAF alignment:"
echo "  Query:  $query_name, length $query_length"
echo "  Target: $target_name, length $target_length"
echo "  Strand: $strand"

if [[ -z "$lw_index_length" ]]; then
    echo "ERROR: Query sequence $query_name was not found in $LW_FAI" >&2
    exit 1
fi

if [[ -z "$sw_index_length" ]]; then
    echo "ERROR: Target sequence $target_name was not found in $SW_FAI" >&2
    exit 1
fi

if [[ "$query_length" != "$lw_index_length" ]]; then
    echo "ERROR: The PAF query does not appear to be the LW assembly." >&2
    echo "PAF query length: $query_length" >&2
    echo "LW index length:  $lw_index_length" >&2
    exit 1
fi

if [[ "$target_length" != "$sw_index_length" ]]; then
    echo "ERROR: The PAF target does not appear to be the SW assembly." >&2
    echo "PAF target length: $target_length" >&2
    echo "SW index length:   $sw_index_length" >&2
    exit 1
fi

echo "PAF orientation verified: query=LW, target=SW."

###############################################################################
# 3. Create the Circos karyotype
###############################################################################

awk 'BEGIN{OFS=" "}
     {
         id="LW_" $1
         print "chr","-",id,id,0,$2,"blue"
     }' "$LW_FAI" > "$KARYOTYPE"

awk 'BEGIN{OFS=" "}
     {
         id="SW_" $1
         print "chr","-",id,id,0,$2,"red"
     }' "$SW_FAI" >> "$KARYOTYPE"

karyotype_count=$(wc -l < "$KARYOTYPE")

echo
echo "Created $KARYOTYPE with $karyotype_count chromosomes/scaffolds."
head -n 5 "$KARYOTYPE"

###############################################################################
# 4. Convert the PAF alignments into Circos links
#
# Important:
# Do not use PAF column 10 divided by column 11 as the identity filter here.
# These alignments contain large gaps. Instead, use Minimap2's dv:f:
# divergence estimate.
###############################################################################

gzip -cd "$PAF" |
awk \
    -v query_prefix="$QUERY_PREFIX" \
    -v target_prefix="$TARGET_PREFIX" \
    -v minimum_span="$MIN_SPAN" \
    -v maximum_divergence="$MAX_DIVERGENCE" \
    -v minimum_mapq="$MIN_MAPQ" '
BEGIN {
    OFS=" "
}
{
    divergence=""
    primary=0

    # Read optional PAF tags.
    for (i=13; i<=NF; i++) {
        if ($i ~ /^dv:f:/) {
            split($i,tag,":")
            divergence=tag[3]+0
        }

        if ($i == "tp:A:P") {
            primary=1
        }
    }

    query_span=$4-$3
    target_span=$9-$8

    # Retain primary, long, high-quality, low-divergence alignments.
    if (primary &&
        divergence != "" &&
        divergence <= maximum_divergence &&
        $12 >= minimum_mapq &&
        query_span >= minimum_span &&
        target_span >= minimum_span) {

        if ($5 == "-") {
            link_color="red_a3"
        } else {
            link_color="blue_a3"
        }

        print query_prefix $1,$3,$4,
              target_prefix $6,$8,$9,
              "color=" link_color
    }
}' > "$LINKS"

link_count=$(wc -l < "$LINKS")

echo
echo "Created $LINKS with $link_count synteny links."
echo "Filters:"
echo "  Minimum span:       $MIN_SPAN bp"
echo "  Maximum divergence: $MAX_DIVERGENCE"
echo "  Minimum MAPQ:       $MIN_MAPQ"
echo "  Primary alignments only"

if [[ "$link_count" -eq 0 ]]; then
    echo "ERROR: No alignments passed the filters." >&2
    exit 1
fi

head -n 5 "$LINKS"

###############################################################################
# 5. Write the Circos configuration
###############################################################################

cat > "$CONFIG" <<'CIRCOS_CONFIG'
<<include etc/colors_fonts_patterns.conf>>

karyotype = karyotype.txt

chromosomes_units           = 1000000
chromosomes_display_default = yes

<ideogram>
radius           = 0.90r
thickness        = 30p
fill             = yes
stroke_color     = black
stroke_thickness = 1p

show_label       = yes
label_radius     = 1.01r
label_size       = 20p
label_parallel   = yes

<spacing>
default = 0.006r
</spacing>
</ideogram>

<ticks>
show_ticks       = yes
show_tick_labels = yes
radius           = 1r
color            = black
thickness        = 2p
multiplier       = 1e-6

<tick>
spacing      = 10u
size         = 8p
show_label   = yes
label_size   = 13p
label_offset = 5p
format       = %d
</tick>
</ticks>

<links>
<link>
file          = links.txt
radius        = 0.85r
bezier_radius = 0.15r
ribbon        = yes
thickness     = 1
</link>
</links>

<image>
dir        = .
file       = Pchalceus_LW_vs_SW.png
png        = yes
svg        = yes
radius     = 1800p
background = white
</image>

<<include etc/housekeeping.conf>>
CIRCOS_CONFIG

echo
echo "Created $CONFIG."

###############################################################################
# 6. Validate chromosome names in the links
###############################################################################

awk '
    NR==FNR {
        valid[$3]=1
        next
    }
    {
        if (!($1 in valid)) {
            print "ERROR: Unknown chromosome in links file: " $1 > "/dev/stderr"
            errors++
        }

        if (!($4 in valid)) {
            print "ERROR: Unknown chromosome in links file: " $4 > "/dev/stderr"
            errors++
        }
    }
    END {
        if (errors)
            exit 1
    }
' "$KARYOTYPE" "$LINKS"

echo "All link chromosome names match the karyotype."

###############################################################################
# 7. Run Circos
###############################################################################

echo
echo "Running Circos..."

circos -conf "$CONFIG"

###############################################################################
# 8. Report the results
###############################################################################

echo
echo "Circos completed successfully."
echo
echo "Output files:"
echo "  ${OUTPUT_PREFIX}.png"
echo "  ${OUTPUT_PREFIX}.svg"

ls -lh \
    "${OUTPUT_PREFIX}.png" \
    "${OUTPUT_PREFIX}.svg"
