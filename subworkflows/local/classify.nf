process qiime2_blast {

    label 'memory_high'
    label 'qiime2'

    publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'
    publishDir "${params.outputdir}/allout", mode: 'copy'

  input:
    path repsep_fasta
    val refseq
    val reftax

  output:
    path 'taxonomy_raw.tsv'

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

    qiime tools export --input-path taxonomy.qza --output-path .
    mv taxonomy.tsv taxonomy_raw.tsv
  """
}

process qiime2_bayes {
    label 'memory_high'
    label 'qiime2'
    publishDir "${params.outputdir}/classify", mode: 'copy'
    publishDir "${params.outputdir}/allout", mode: 'copy'

  input:
    path repsep_fasta
    path model

  output:
    path 'taxonomy_raw.tsv', emit: taxonomy_raw

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
  mv taxonomy.tsv taxonomy_raw.tsv
  """
}

process clean_taxonomy_tsv {
  publishDir "${params.outputdir}/allout", mode: 'copy'

  container "amancevice/pandas:1.3.5"

  input:
    path taxa_raw
  
  output:
    path 'taxonomy.tsv', emit: taxtab
  """
    # ${params.dada.chimera_alg} 
    parse_qiime.py ${taxa_raw} taxonomy.tsv qiime2_silva
  """
}

workflow classify_reads {
  take: 
    fasta_reads
    bayesmodel
    refseq
    reftax
// Last two for BLAST search
  main:
    if (bayesmodel != null) {
      qiime2_bayes(fasta_reads, bayesmodel)
      clean_taxonomy_tsv(qiime2_bayes.out.taxonomy_raw)
      cout = clean_taxonomy_tsv.out.taxtab
    }
    else if (refseq != null){
      qiime2_blast(refseq, reftax)
      cout = qiime2_blast.out
    } else {

    }
  emit:
    taxtab = cout
}

// vi: ft=groovy
