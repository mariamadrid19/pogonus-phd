#!/usr/bin/env python3
from Bio.Blast import NCBIXML

E_VALUE_THRESH = 0.04

with open("augustus.blastout.swiss.xml") as blast_file, open("augustus_parsed_blastout.tsv", "w") as output_file:
    for record in NCBIXML.parse(blast_file):
        if record.alignments:
            for j, align in enumerate(record.alignments):
                for i, hsp in enumerate(align.hsps):
                    if (j == 0) & (i == 0):
                        output_file.write(f"{j}\t{i}\t{record.query}\t{hsp.expect}\t{align.title}\n")

with open("braker.blastout.swiss.xml") as blast_file, open("braker_parsed_blastout.tsv", "w") as output_file:
    for record in NCBIXML.parse(blast_file):
        if record.alignments:
            for j, align in enumerate(record.alignments):
                for i, hsp in enumerate(align.hsps):
                    if (j == 0) & (i == 0):
                        output_file.write(f"{j}\t{i}\t{record.query}\t{hsp.expect}\t{align.title}\n")
