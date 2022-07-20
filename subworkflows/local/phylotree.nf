/*
  Building phylogenetic tree
*/
process qiime2_roottree {

  label 'big_mem'
  label 'qiime2'
  publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'
  publishDir "${params.outputdir}/allout", mode: 'copy'
  errorStrategy "ignore"

  input:
    path repsep_fasta
    path asvtab
    val roottree
  output:
    path 'rooted-tree.qza'
    path 'rooted-tree.nwk'
    path 'filtered_table.tsv'
    path 'removed_table.biom'

  shell:
  """
  export MPLCONFIGDIR=/tmp/abc123
  qiime tools import --input-path ${repsep_fasta} \
    --output-path sequences.qza \
    --type 'FeatureData[Sequence]'

  biom convert -i ${asvtab} -o asv.biom --to-hdf5

  qiime tools import \
    --input-path asv.biom \
    --type 'FeatureTable[Frequency]' \
    --output-path table.qza

  qiime fragment-insertion sepp \
    --i-representative-sequences ./sequences.qza \
    --i-reference-database ${roottree} \
    --o-tree ./rooted-tree.qza \
    --o-placements ./tree_placements.qza \
    --p-threads ${task.cpus}

  qiime fragment-insertion filter-features \
    --i-table table.qza \
    --i-tree ./rooted-tree.qza \
    --o-filtered-table filtered_table.qza \
    --o-removed-table removed_table.qza
  
  qiime tools export --input-path filtered_table.qza --output-path .
  biom convert -i feature-table.biom -o filtered_table.tsv --to-tsv

  rm feature-table.biom

  qiime tools export --input-path removed_table.qza --output-path .
  mv feature-table.biom removed_table.biom

  qiime tools export --input-path rooted-tree.qza --output-path  .
  mv tree.nwk rooted-tree.nwk
  """
}

process qiime2_roottree_mafft {
    label 'med_mem'
    label 'qiime2'
    publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'
    publishDir "${params.outputdir}/allout", mode: 'copy'
  input:
    path repsep_fasta
  output:
    path 'rooted-tree.qza'
    path 'rooted-tree.nwk'

  shell:
  """
  export MPLCONFIGDIR=/tmp/xyz
  qiime tools import --input-path ${repsep_fasta} \
    --output-path sequences.qza \
    --type 'FeatureData[Sequence]'

  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences sequences.qza \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza \
    --p-n-threads ${task.cpus} || true

  qiime tools export --input-path rooted-tree.qza --output-path  .
  mv tree.nwk rooted-tree.nwk
  """
}


workflow qiime2_tree {
  take:
    fasta
    model

  main:
    if (params.root-mafft){
    }
    else {
    }

  emit:
    qiime2_roottree_mafft.out
}
