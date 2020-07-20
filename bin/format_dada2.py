#!/usr/bin/env python

""" Format DADA output.

This script extract the row-name (FASTA-sequence) of DADA2 output and move it into another files.
Also change the sample name (since DADA2 use first read as a filename)
"""

import sys
import re
import os
import hashlib
import pandas as pd

INPUT=sys.argv[1]
# If pass the manifest, rename columns as well

try:
    MANIFEST=sys.argv[2]
    with open(MANIFEST) as fh:
        mandf = pd.read_csv(MANIFEST, sep="\t")
        # Validate and make sure that there is no duplicated data.
        # First, forward and reverse should map to the same file.
        renmapper = dict(zip(mandf["filename"], mandf["sampleID"]))
except IndexError:
    MANIFEST=None


# Reformat table using MD5sum and output sequence for later use.
a = pd.read_csv(INPUT, sep="\t", index_col=0)
if MANIFEST:
    a = a.rename(mapper=renmapper, axis=1)
ASV_SEQ = a.index
indexname = a.index.name
id_md5 = [hashlib.md5(i.encode('utf-8')).hexdigest() for i in ASV_SEQ]

# Create fasta file
with open("repsep.fasta", "w") as fho:
    for header, fasta in zip(id_md5, ASV_SEQ):
        fho.write(">{}".format(header))
        fho.write(os.linesep)
        fho.write(fasta)
        fho.write(os.linesep)

# Create table file with md5sum
a.index = id_md5
a.index.name = indexname
# Rename
if MANIFEST:
    a = a.rename(columns=renmapper)

with open("asv.tab", "w") as fho:
    a.to_csv(fho, sep="\t")
