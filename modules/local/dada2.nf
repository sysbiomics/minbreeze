// DADA2 for single-read / merged data
process dada2_single {

  label 'process_medium'

  conda (params.enable_conda ? "${projectDir}/envs/dada2.yaml" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' :
      'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

  publishDir "${params.outputdir}/dada2", mode: 'copy'
  publishDir "${params.outputdir}/allout", mode: 'copy'

  input:
    path(reads)

  output:
    path 'dadaraw.tsv', emit: dada_raw_table
    path 'track.tsv', emit: track
    path 'qualplot.pdf'

  """
    mkdir -p merged
    mv ${reads} merged
    run_dada_single.R merged dadaraw.tsv track.tsv mm 0 0 2.0 2 500 ${params.dada.chimera_alg} ${params.dada.chimera_fol} ${task.cpus} 1000000 ${params.dada.pool} NULL 16
  """
}

process dada2_pair {

  label 'process_medium'

  conda (params.enable_conda ? "${projectDir}/envs/dada2.yaml" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' :
      'quay.io/biocontainers/bioconductor-dada2:1.22.0--r41h399db7b_0' }"

  publishDir "${params.outputdir}/dada2", mode: 'copy'
  publishDir "${params.outputdir}/allout", mode: 'copy'

  input:
    path freads
    path rreads

  output:
    path 'dadaraw.tsv', emit: dada_raw_table
    path 'track.tsv', emit: track
    path 'qualplotF.pdf'
    path 'qualplotR.pdf'

  """
  mkdir -p fwd
  mkdir -p rev
  mv ${freads} fwd
  mv ${rreads} rev
  run_dada_paired.R fwd rev dadaraw.tsv track.tsv ff rf 0 0 0 0 2.0 2.0 2 ${params.dada.chimera_alg} ${params.dada.chimera_fol} ${task.cpus} 1000000 ${params.dada.pool} ${params.minOverlap}
  """
}


process export_dada2tsv {

  publishDir "${params.outputdir}/dada2", mode: 'copy'
  publishDir "${params.outputdir}/allout", mode: 'copy'

  conda (params.enable_conda ? "${projectDir}/envs/pandas.yaml" : null)
  container "amancevice/pandas:1.3.5"

  input:
    path dada_raw

  output:
    path 'asv.tab', emit: asvtab
    path 'repsep.fasta', emit: repsep
    val "NOT QIIME2", emit: source

  """
    format_dada2.py ${dada_raw}
  """
}
