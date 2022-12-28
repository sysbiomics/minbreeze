TODO
====

1. Convert taxonomy into .tsv.
2. Taxonomy classification
  2.1 Check how naive-bayes match since I am not sure if there is a identity limit
  2.2 Implement BLAST and Vsearch as an alternative way
3. Run minpath
4. Cache all model from QIIME2 instead of redownloading everytime.
5. Report read with each steps


TO_FIX
=====


## TODO
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
