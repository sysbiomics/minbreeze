#!/usr/bin/python

import sys
import re
import os
import hashlib
import pandas as pd


INPUT=sys.argv[1]
# Reformat table using MD5sum and output sequence for later use.
a = pd.read_csv(INPUT, sep="\t")
ASV_SEQ = a["#OTU ID"]
    
id_md5 = [hashlib.md5(i.encode('utf-8')).hexdigest() for i in ASV_SEQ]

# Create fasta file
with open("repsep.fasta", "w") as fho:
    for header, fasta in zip(id_md5, ASV_SEQ):
        fho.write(">{}".format(header))
        fho.write(os.linesep)
        fho.write(fasta)
        fho.write(os.linesep)

# Create table file with md5sum
a["#OTU ID"] = id_md5
# Remove "_R?.fastq.gz from name"

with open("asv.tab", "w") as fho:
    a.to_csv(fho, sep="\t", index=False)
