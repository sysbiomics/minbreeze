/*
*  DAF pipeline
*
*/


def helpMessage() {
  log.info"""
  --inputdir     Path to input data
  --outputdir    Path to output data
  """.stripIndent()
}


println("${params.getClass()}")
// Show help message
if (params.help){
    helpMessage()
    exit 0
}

pmatcher = "*_R{1,2}.fastq.gz"
matcher = params.inputdir + "/" + pmatcher

println("${matcher}" )
params.fwdprimer = "CCTACGGGNGGCWGCAG"
params.revprimer = "GACTACHVGGGTATCTAATCC"

/*
*  Create a channel for input read files
*/

Channel
  .fromFilePairs(matcher)
  .ifEmpty { exit 1, "Cannot find anything here chief"}
  .into { ch_read_pairs; ch_read_pairs_fastqc}

/*
  Start with FASTQC to report initial result
*/
process fastqc {
  publishDir "${params.outputdir}/fastqc", mode: 'copy'
  input:
  set val(pair_id), file(reads) from ch_read_pairs_fastqc
  output:
  file "*.html"
  file "*.zip"

  shell:
  """
  fastqc ${reads}
  """
}

/*
process multiqc {
}
*/

/*
  Trimming adapter and quality with BBDUK

  Using DADA2 is a bit hard to automate, so we need to compromise on this
*/
literals=params.fwdprimer + "," + params.revprimer

process bbduk {

  publishDir "${params.outputdir}", mode: 'copy'
  input:
  set val(pair_id), file(reads) from ch_read_pairs

  output:
  file("trimmed/${reads[0]}") into ch_dada2forReads
  file("trimmed/${reads[1]}") into ch_dada2revReads
  file("trimmed/*.log")

  """
  mkdir -p trimmed
  echo ${reads[0]}
  bbduk.sh -Xmx1g in=${reads[0]} in2=${reads[1]} out1=trimmed/${reads[0]} \
  out2=trimmed/${reads[1]} \
  literal=${literals} \
  ktrim=l k=13 mink=6 hdist=1 qtrim=r trimq=15 \
  minlength=150 restrictleft=25 copyundefined=T stats=trimmed/stat_${pair_id}.log \
  2> run_${pair_id}.log
  """
}


process dada2 {
  label 'big_cpu'

  publishDir "${params.outputdir}/dada2", mode: "copy"
  input:
    file freads from ch_dada2forReads.collect()
    file rreads from ch_dada2revReads.collect()

  output:
    file("dadaraw.tsv") into ch_dada2ext
    file("track.tsv")
    file("asv.tab")
    file("repsep.fasta")
  """
  mkdir -p fwd
  mkdir -p rev
  mv ${freads} fwd
  mv ${rreads} rev
  run_dada_paired.R fwd rev dadaraw.tsv track.tsv ff rf 0 0 0 0 2.0 2.0 2 consensus 1.0 ${task.cpus} 250000
  """
}

/*
process dada2extract {
  publishDir "${params.outputdir}/dada2", mode: "copy"

  input:
    file rawasv from ch_dada2ext

  output:
    file("asv.tab")
    file("repsep.fasta")

  """
  format_dada2.py ${rawasv}  # Produce asv.tab and repsep.fasta
  """
}
*/

