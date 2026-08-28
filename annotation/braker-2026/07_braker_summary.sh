#!/bin/bash -l
#SBATCH --cluster=wice
#SBATCH --job-name=braker_summary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=12:00:00
#SBATCH --output=braker_summary.%j.out
#SBATCH --account=lp_svbelleghem

DIR=./annotation/braker
GTF="$DIR/braker.gtf"
AA="$DIR/braker.aa"
CDS="$DIR/braker.codingseq"

[[ -s "$GTF" ]] || {
    echo "ERROR: $GTF is missing or empty" >&2
    exit 1
}

echo "BRAKER annotation summary"
echo "Directory: $DIR"
echo

awk -F '\t' '
function attribute(name, text,    pattern,value) {
    pattern = name "[ =\"]+[^;\" ]+"
    if (match(text, pattern)) {
        value = substr(text, RSTART, RLENGTH)
        sub("^" name "[ =\"]+", "", value)
        return value
    }
    return ""
}

/^#/ { next }

NF == 9 {
    feature[$3]++
    scaffold[$1] = 1

    gene = attribute("gene_id", $9)
    transcript = attribute("transcript_id", $9)

    if (gene != "")
        genes[gene] = 1

    if (transcript != "") {
        transcripts[transcript] = 1

        if (gene != "") {
            pair = gene SUBSEP transcript
            gene_transcript[pair] = 1
        }
    }

    if ($3 == "CDS" && transcript != "")
        transcript_has_cds[transcript] = 1
}

END {
    for (x in genes)
        ngenes++

    for (x in transcripts)
        ntranscripts++

    for (x in scaffold)
        nscaffolds++

    for (x in gene_transcript) {
        split(x, p, SUBSEP)
        isoforms[p[1]]++
    }

    for (x in transcript_has_cds)
        coding_transcripts++

    single_isoform = 0
    multi_isoform = 0
    maximum_isoforms = 0

    for (g in genes) {
        n = isoforms[g] + 0
        if (n <= 1)
            single_isoform++
        else
            multi_isoform++

        if (n > maximum_isoforms)
            maximum_isoforms = n
    }

    print "Genes:", ngenes + 0
    print "Transcripts:", ntranscripts + 0
    print "Coding transcripts:", coding_transcripts + 0
    print "Mean transcripts per gene:", \
        ngenes ? ntranscripts / ngenes : 0
    print "Single-isoform genes:", single_isoform
    print "Multi-isoform genes:", multi_isoform
    print "Maximum isoforms for one gene:", maximum_isoforms
    print "Annotated scaffolds/contigs:", nscaffolds
    print ""
    print "Feature rows:"
    print "  genes:", feature["gene"] + 0
    print "  transcripts/mRNAs:", \
        feature["transcript"] + feature["mRNA"]
    print "  exons:", feature["exon"] + 0
    print "  CDS segments:", feature["CDS"] + 0
    print "  start codons:", feature["start_codon"] + 0
    print "  stop codons:", feature["stop_codon"] + 0
}' "$GTF"

if [[ -s "$AA" ]]; then
    echo
    echo "Predicted protein statistics:"

    awk '
    function finish() {
        if (!active)
            return

        n++
        lengths[n] = length(sequence)
        total += lengths[n]

        if (lengths[n] < 100)
            under100++
        if (lengths[n] > maximum)
            maximum = lengths[n]

        sequence = ""
    }

    /^>/ {
        finish()
        active = 1
        next
    }

    {
        line = $0
        gsub(/[[:space:]]/, "", line)
        gsub(/\*$/, "", line)
        sequence = sequence line
    }

    END {
        finish()

        for (i = 1; i <= n; i++)
            for (j = i + 1; j <= n; j++)
                if (lengths[i] > lengths[j]) {
                    temp = lengths[i]
                    lengths[i] = lengths[j]
                    lengths[j] = temp
                }

        if (n % 2)
            median = lengths[(n + 1) / 2]
        else if (n)
            median = (lengths[n / 2] + lengths[n / 2 + 1]) / 2

        print "  protein sequences:", n
        print "  total amino acids:", total
        printf "  mean protein length: %.1f aa\n", n ? total/n : 0
        printf "  median protein length: %.1f aa\n", median
        print "  longest protein:", maximum, "aa"
        print "  proteins shorter than 100 aa:", under100 + 0
    }' "$AA"
else
    echo
    echo "WARNING: $AA is absent or empty; protein statistics skipped."
fi
