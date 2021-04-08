# MinBreeze

The analysis pipeline use in for pair-end 16s dataset.

## Why this and not QIIME2?
QIIME2 does an excellent job in making everything into the pipeline (trackable process, parallel run), but the effort making everything accessable end up hidden all the process. Good for the beginner I guess, but it becomes a bit annoying when you want to tinker

## What does this do then?
16s ofcourse

## How to run the pipeline


## TODO
1. Report read with each steps
2. Taxonomy classification
  2.1 Check how naive-bayes match since I am not sure if there is a identity limit
  2.2 Implement BLAST and Vsearch as an alternative way
3. Phylogenetic tree
  4.1 https://github.com/qiime2/q2-fragment-insertion/blob/master/q2_fragment_insertion/_insertion.py
4. Validation of result
  4.1 Internal checking with multiqc. (https://multiqc.info/docs/#custom-content)
  4.2 Final result checking with external data (avaialble or insilico generate)
  

# Reference
Add protocal https://help.ezbiocloud.net/16s-mtp-protocol-for-illumina-iseq-100/
How to beta diversity with https://astrobiomike.github.io/amplicon/dada2_workflow_ex
How to make q2-fragment insertion https://forum.qiime2.org/t/q2-fragment-insertion-reference-info/4360 , with reference here https://raw.githubusercontent.com/qiime2/q2-fragment-insertion/master/bin/import_references.py
and the SEPP https://github.com/qiime2/q2-fragment-insertion/blob/master/q2_fragment_insertion/_insertion.py
