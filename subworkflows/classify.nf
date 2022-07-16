process qiime2_blast {
    label 'big_mem'
    label 'qiime2'
    publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'
    publishDir "${params.outputdir}/allout", mode: 'copy'

  input:
    path repsep_fasta
    val refseq
    val reftax

  output:
    path 'taxonomy.tsv'

    """

  qiime tools import --input-path ${repsep_fasta} \
    --output-path sequences.qza \
    --type 'FeatureData[Sequence]'

  qiime feature-classifier classify-consensus-blast \
    --i-query sequences.qza \
    --i-reference-reads ${refseq} \
    --i-reference-taxonomy ${reftax} \
    --o-classification taxonomy.qza \
    --p-strand 'plus' \
    --p-unassignable-label unassigned

  qiime tools export  --input-path taxonomy.qza --output-path .
  """
}

process qiime2_bayes {
    label 'big_mem'
    label 'qiime2'
    publishDir "${params.outputdir}/classify", mode: 'copy'
    publishDir "${params.outputdir}/allout", mode: 'copy'

  input:
    path repsep_fasta
    path model

  output:
    path 'taxonomy.tsv'

    """
  qiime tools import --input-path ${repsep_fasta} \
    --output-path sequences.qza \
    --type 'FeatureData[Sequence]'

  qiime feature-classifier classify-sklearn \
    --i-classifier ${model} \
    --p-n-jobs ${task.cpus} \
    --i-reads  sequences.qza \
    --o-classification taxonomy.qza

  qiime tools export --input-path taxonomy.qza --output-path .
  """
}

// The cleanest way would be using switch.
workflow classify_reads {
  take: 
    fasta_reads
    models
  main:
    switch (params.sklearn) {
        case true:
                qiime2bayes()
                break
        default:
                qiime2blast()
                break
        }
  emit:
    qiime2bayes.out
}


// vi: ft=groovy
