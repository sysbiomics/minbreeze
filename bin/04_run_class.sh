#!/bin/bash

set -e
#conda activate qiime2-2019.7

source activate qiime2-2019.1
cd $PBS_O_WORKDIR

# qiime tools import --input-path dadaoutput/repsep.fasta  --output-path sequences.qza   --type 'FeatureData[Sequence]'

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences sequences.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime feature-classifier classify-sklearn \
  --i-classifier ../../res/gg99_341f-805r_classifier.qza \
  --p-n-jobs 15 \
  --i-reads  sequences.qza \
  --o-classification taxonomy_gg.qza
