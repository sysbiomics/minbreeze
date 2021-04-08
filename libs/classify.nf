process qiime2_blast {
    label 'mod_cpu'
    publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'

  input:
    path repsep_fasta
    val refseq
    val reftax

  output:
    path 'taxonomy.qza'
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

  extract_fl () {
    local zipfile=\$1
    local lookfor=\$2
    local location_z=\$(unzip -l \$zipfile  | gawk '{print \$4}' | grep \$lookfor)  # Hack, but work
    echo \$location_z
    # Unzip extract a lot of file, so force overwrite is neccessary
    unzip -o -j \$zipfile \$location_z
  }

  extract_fl taxonomy.qza taxonomy.tsv
  """
}

process qiime2_bayes {
    label 'mod_cpu'
    publishDir "${params.outputdir}/classify", mode: 'copy'

  input:
    path repsep_fasta
    val model

  output:
    path 'taxonomy.qza'
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

  extract_fl () {
    local zipfile=\$1
    local lookfor=\$2
    local location_z=\$(unzip -l \$zipfile  | gawk '{print \$4}' | grep \$lookfor)  # Hack, but work
    echo \$location_z
    # Unzip extract a lot of file, so force overwrite is neccessary
    unzip -o -j \$zipfile \$location_z
  }

  extract_fl taxonomy.qza taxonomy.tsv
  """
}

// The cleanest way would be using switch.
workflow classify_reads {
  take: fasta_reads
  main:
    switch (params.qtrim) {
    case true:
            dettrim(chan, params.fwdprimerlen, params.revprimerlen)
            dada2(dettrim.out.fwd_reads.collect(), dettrim.out.rev_reads.collect(), params.fwdlen, params.revlen, params.fwdprimerlen, params.revprimerlen)
            break
    default:
        dada2(something.fwd.collect(), something.rev.collect(), params.fwdlen, params.revlen, params.fwdprimerlen, params.revprimerlen)
            break
    }
  emit:
    classified_table
}

// vi: ft=groovy
