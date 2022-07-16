// Multiple way to clean reads

process fastqc {
    publishDir "${params.outputdir}/fastqc", mode: 'copy'

  input:
    tuple val(pair_id), file(reads)
  output:
    file '*.html'
    file '*.zip'

  shell:
    """
  fastqc ${reads[0]}
  fastqc ${reads[1]}
  """
}

process QCTrim {
    """
    bbduk.sh -Xmx1g in=t_trimleft/${reads[0]} in2=t_trimleft/${reads[1]} \
    out1=trimmed/${reads[0]} out2=trimmed/${reads[1]} \
    qtrim=r trimq=15 \
    minlength=150 stats=trimmed/stat_${pair_id}.log \
    2> trimmed/run_${pair_id}.log
    echo "Deterministic trim with ${params.fwdprimerlen} ${params.revprimerlen}" > trimmed/seqtk.log
    """
}

workflow QCTrim {
  take: reads
  main:
    QCTrim(reads)
  emit:
    QCTrim.out
}
// vi: ft=groovy
