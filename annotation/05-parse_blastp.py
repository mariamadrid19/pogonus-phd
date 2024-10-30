#!/usr/bin/env python3
from Bio.Blast import NCBIXML

E_VALUE_THRESH = 0.04

with open("blastout.swiss.h.xml") as blast_file, open("parsed_blastout.tsv", "w") as output_file:
    for record in NCBIXML.parse(blast_file):
        if record.alignments:
            for j, align in enumerate(record.alignments):
                for i, hsp in enumerate(align.hsps):
                    if (j == 0) & (i == 0):
                       
