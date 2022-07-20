TODO
====

1. Convert taxonomy into .tsv.
2. Run PICRUSt2
3. Run minpath


TOFIX
=====
1. Singularity has problem with symlink, which mean the nextflow will have problem with it.
  - Point, MAFFT cannot function at all.


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
