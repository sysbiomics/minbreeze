/*
 *  DAF pipeline
 *
 */

nextflow.preview.dsl = 2
ANSI_RED = "\033[1;31m"
ANSI_GREEN = "\033[1;32m"
ANSI_BLUE = "\033[1;34m"
ANSI_RESET = "\033[0m"

/*
 * Default pipeline parameters
 */
import java.nio.file.Files
import java.nio.file.Paths
params.help = ''

nextflowVersion = '>=20.04.0'        // 1.2 or later

def makeabs(filename) {
    return new File(filename).getCanonicalPath()
}

def required(params) {
    // Error if parameters is not available
    throw new Exception('Param not exists')
}

/*
Redefine and clean parameters
*/
// get counts of found fastq files
print params
reads = "${params.inputdir}/${params.fqpattern}"
readcounts = file(reads).size()

params.qtrim=true
params.modeltax = new File(params.modeltax).getCanonicalPath()
params.manifest = new File(params.manifest).getCanonicalPath()
params.refseqabs = new File(params.refseq).getCanonicalPath()
params.reftaxabs = new File(params.reftax).getCanonicalPath()
params.roottree = new File(params.modeltree).getCanonicalPath() // Root tree for sepp

/*
switch(params.classifier){
  case "BLAST":
    include { qiime2_blast as classifier } from './libs/classify.nf'
    def x = 20
  case "QIIME":
    include { qiime2_bayes_ as classifier } from './libs/classify.nf'
  default:
    exit(1)
}
*/

/* End of pipeline parameter
*/

log.info """
  =================================
  NXF-minbreeze
  =================================

  Parameters:
  -------------------------------------------
  --inputdir        : ${params.inputdir}
  --fqpattern       : ${params.fqpattern}

  Runtime data:
  -------------------------------------------
  Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
  Container:              ${ANSI_GREEN}${workflow.container}${ANSI_RESET}
  Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
  Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
  Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
  Fastq files:            ${ANSI_GREEN}${ readcounts } files found${ANSI_RESET}
  """
         .stripIndent()


def helpMessage() {
    log.info """
  =================================
  SALSA
  =================================

  Usage:
  ---------------------------------

  --inputdir     Path to input data
  --outputdir    Path to output data
  --config       Path to config
  """.stripIndent()
}

// Help message
if (params.help) {
    helpMessage()
    exit 0
}

process qc_fastqc {
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

/*
process trunc_seqtk {
}
*/

/*
 * Trimming primer and quality. We don't use DADA2 for this because we might
 * want to try other trimming solution instead.
*/
process trim_seqtkbbduk {
    publishDir "${params.outputdir}/trimmed", mode: 'link'

  input:
    tuple val(pair_id), file(reads)
    val(truncleft1)
    val(truncleft2)
//    val(qtrim)
  output:
    path "trimmed/${reads[0]}", emit: fwd_reads
    path "trimmed/${reads[1]}", emit: rev_reads
    path 'trimmed/*.log'

//  script:
//    bbduk_qtrim= qtrim ? "qtrim=r trimq=$qtrim" : ""

  shell:
  """
	mkdir -p t_trunleft
	mkdir -p trimmed
	echo ${pair_id}

	bbduk.sh ftl=${truncleft1} ${reads[0]} in=${reads[0]} out=t_trunleft/${reads[0]}
	bbduk.sh ftl=${truncleft2} ${reads[1]} in=${reads[1]} out=t_trunleft/${reads[1]}

	bbduk.sh -Xmx1g in=t_trunleft/${reads[0]} in2=t_trunleft/${reads[1]} \
	  out1=trimmed/${reads[0]} out2=trimmed/${reads[1]} \
	  qtrim=r trimq=15 \
	  minlength=150 stats=trimmed/stat_${pair_id}.log \
	  2> trimmed/run_${pair_id}.log
	  echo "Deterministic trim with ${params.fwdprimerlen} ${params.revprimerlen}" > trimmed/seqtk.log
  """
}

/* Denoise
*/
process asvcall_dada2 {
    label 'big_cpu'

    publishDir "${params.outputdir}/dada2", mode: 'copy'
  input:
    path freads
    path rreads
    val fwdlen
    val revlen
    val trimleft
    val trimright

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
  run_dada_paired.R fwd rev dadaraw.tsv track.tsv ff rf ${fwdlen} ${revlen} ${trimleft} ${trimright} 2.0 2.0 2 ${params.dada.chimera_alg} ${params.dada.chimera_fol} ${task.cpus} 1000000 ${params.dada.pool}
  """
}

process export_dada2tsv {
    publishDir "${params.outputdir}/dada2", mode: 'copy'

  input:
    path rawasv

  output:
    path 'asv.tab'
    path 'repsep.fasta', emit: repsep
  script:
    def manifest = params.manifestabs ? params.manifestabs : ''
    """
    format_dada2.py ${rawasv} ${manifest}
    """
}

/*
  Building phylogenetic tree
*/
process qiime2_roottree {
    label 'big_cpu'
    publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'
  input:
    path repsep_fasta
    val roottree
  output:
    path 'rooted-tree.qza'
    path 'rooted-tree.nwk'

  shell:
  """
  qiime tools import --input-path ${repsep_fasta} \
    --output-path sequences.qza \
    --type 'FeatureData[Sequence]'

  #qiime phylogeny align-to-tree-mafft-fasttree \
  #  --i-sequences sequences.qza \
  #  --o-alignment aligned-rep-seqs.qza \
  #  --o-masked-alignment masked-aligned-rep-seqs.qza \
  #  --o-tree unrooted-tree.qza \
  #  --o-rooted-tree rooted-tree.qza

  qiime fragment-insertion sepp \
    --i-representative-sequences ./sequences.qza \
    --i-reference-database ${roottree} \
    --o-tree ./rooted-tree.qza \
    --o-placements ./tree_placements.qza \
    --p-threads ${task.cpus}  # update to a higher number if you can

  extract_fl () {
    local zipfile=\$1
    local lookfor=\$2
    local location_z=\$(unzip -l \$zipfile  | gawk '{print \$4}' | grep \$lookfor)  # Hack, but work
    echo \$location_z
    # Unzip extract a lot of file, so force overwrite is neccessary
    unzip -o -j \$zipfile \$location_z
  }

  extract_fl rooted-tree.qza tree.nwk
  mv tree.nwk rooted-tree.nwk
  """
}

include { qiime2_bayes as classifier } from './libs/classify.nf'

/*
*  Final workflow
*/

workflow {
    chan = Channel
    .fromFilePairs(reads)
    .ifEmpty { exit 1, 'No input files' }
    qc_fastqc(chan)
    if (params.qtrim == true) {
        trim_seqtkbbduk(chan, params.fwdprimerlen, params.revprimerlen)
        asvcall_dada2(trim_seqtkbbduk.out.fwd_reads.collect(), trim_seqtkbbduk.out.rev_reads.collect(), 0, 0, 0, 0)}
    else {
        chan.multiMap { it ->
            fwd: it[1][0]
            rev: it[1][1]}.set { rawseqs }
        asvcall_dada2(rawseqs.fwd.collect(), rawseqs.rev.collect(), params.fwdlen, params.revlen, params.fwdprimerlen, params.revprimerlen)
        }

    export_dada2tsv(asvcall_dada2.out.dada_raw_table)
    qiime2_roottree(export_dada2tsv.out.repsep, params.roottree)
    classifier(export_dada2tsv.out.repsep, params.modeltax)
    }

// Export options use in this into json.
import groovy.json.JsonBuilder

workflow.onComplete {
    // serialize
    def builder = new JsonBuilder()
    builder(params)
    def output_f = new File("${params.outputdir}/pipeline.txt")
    output_f.withWriter { w -> w << builder.toString() }
}

// vi: ft=groovy
