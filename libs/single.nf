/*
  Single version
*/
process fastqc_single {
  input:
    file(read)
  output:
    file "*.html"
    file "*.zip"

  shell:
  """
  fastqc ${read}
  """
}

/*  Determistic trimming, with single
*/
process dettrim_single {

  input:
    file(read)

  output:
    path "trimmed/${read}", emit: read
    path "trimmed/*.log"

  """
  mkdir -p t_trimleft
  mkdir -p trimmed
  echo ${read}
  seqtk trimfq -b ${params.fwdprimerlen} ${read} | gzip > t_trimleft/${read}

  bbduk.sh -Xmx1g in=t_trimleft/${read} \
  out=trimmed/${read} \
  qtrim=r trimq=15 \
  minlength=150 stats=trimmed/stat_${read}.log \
  2> trimmed/run_${read}.log
  echo "Deterministic trim with ${params.fwdprimerlen}" > trimmed/seqtk.log
  """
}
