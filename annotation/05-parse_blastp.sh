#!/bin/bash -l
#SBATCH --cluster=genius
#SBATCH --job-name parse_blastP
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --time=48:00:00
#SBATCH -o parse_blastP.%j.out
#SBATCH -A lp_svbelleghem

module load Biopython/1.83-foss-2023a
python3
from Bio.Blast import NCBIXML
 
E_VALUE_THRESH = 0.04
 
for record in NCBIXML.parse(open("blastout.swiss.h.xml")):
        if record.alignments:
                #print("query: %s" % record.query[:100])
                for j,align in enumerate(record.alignments):
                        for i,hsp in enumerate(align.hsps):
                                        if (j==0) & (i==0):
                                                print("%s\t%s\t%s\t%s\t%s" % (j,i,record.query,hsp.expect,align.title))
