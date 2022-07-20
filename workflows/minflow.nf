/*
 *  Minflow pipeline
 *
 */

nextflow.enable.dsl=2

ANSI_RED = "\033[1;31m"
ANSI_GREEN = "\033[1;32m"
ANSI_BLUE = "\033[1;34m"
ANSI_RESET = "\033[0m"

/*
 * Default pipeline parameters
 */
params.help = ''

/*
Redefine and clean parameters
*/
// get counts of found fastq files
// print params
// reads = "${params.inputdir}/${params.fqpattern}"
// readcounts = file(reads).size()

// Tree param
params.rootmafft = true
params.roottree = false

/* End of pipeline parameter
*/

def helpMessage() {
    log.info """
  =================================
  NXF-MINBREEZE
  =================================

  Usage:
  ---------------------------------

  --input     Path to input data
  --outputdir    Path to output data
  --config       Path to config
  """.stripIndent()
}

// Help message
if (params.help) {
    helpMessage()
    exit 0
}

// Report
log.info """
  =================================
  NXF-minbreeze
  =================================

  Parameters:
  -------------------------------------------
  --input        : ${params.input}
  --minOverlap      : ${params.minOverlap}
  --trimfront1      : ${params.trimfront1}
  --trimfront2      : ${params.trimfront2}
  --trimtail1       : ${params.trimtail1}
  --trimtail2       : ${params.trimtail2}
  --modeltax        : ${params.modeltax}
  --modeltree       : ${params.modeltree}

  Runtime data:
  -------------------------------------------
  Running with profile:   ${ANSI_GREEN}${workflow.profile}${ANSI_RESET}
  Container:              ${ANSI_GREEN}${workflow.container}${ANSI_RESET}
  Running as user:        ${ANSI_GREEN}${workflow.userName}${ANSI_RESET}
  Launch dir:             ${ANSI_GREEN}${workflow.launchDir}${ANSI_RESET}
  Base dir:               ${ANSI_GREEN}${baseDir}${ANSI_RESET}
  """.stripIndent()

//
//
//

include { picrust2                               } from '../modules/local/picrust'

//
// Workflow
//

include { INPUT_CHECK                            } from '../subworkflows/local/input_check'
include { qiime2_bayes as classifier             } from '../subworkflows/local/classify'
include { qiime2_roottree_mafft; qiime2_roottree } from '../subworkflows/local/phylotree'


/*
 * Trimming primer and quality. We don't use DADA2 for this because user 
 * might want to try other trimming solution instead.
*/
process fastp {

    publishDir "${params.outputdir}/fastp", mode: 'link'

  input:
    tuple val(meta), path(reads)
    val(truncleft1)
    val(truncleft2)
    val(truncright1)
    val(truncright2)
  output:
    tuple val(meta), path('*.trim.fastq.gz'), emit: reads
    path 'logs/*'

  script:
    def qtrim_tail = ((params.trimtail1 == 0) && (params.trimtail2 == 0)) ? "-3" : ""
  """
    mkdir logs
    fastp -i ${reads[0]} -I ${reads[1]} \
      -f ${truncleft1} -F ${truncleft2} -t ${truncright1} -T ${truncright2} \
      -A -q 15 -l 210 ${qtrim_tail} \
      -o ${meta.id}_R1.trim.fastq.gz -O ${meta.id}_R2.trim.fastq.gz \
      -j logs/${meta.id}.json -h logs/${meta.id}.html
  """
}

/* Merge read
*/
process flash {
  publishDir "${params.outputdir}/flash", mode: 'link'

  input:
    tuple val(meta), path(reads)
  output:
    tuple val(meta), path("${meta.id}.fastq.gz"), emit: reads
    path('logs/*.log'), emit: read

  """
    mkdir -p logs
    flash -m ${params.minOverlap} -M 70 ${reads[0]} ${reads[1]} -o ${meta.id} -t ${task.cpus} -z 2>&1 | tee logs/${meta.id}.log
    # Remove .extendedFrags.fastq.gz part
    mv ${meta.id}.extendedFrags.fastq.gz ${meta.id}.fastq.gz

  """
}

/* Denoise
*/
process dada2 {

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
    run_dada_single.R merged dadaraw.tsv track.tsv mm 0 0 2.0 2 500 consensus 1.0 8 1000000 pseudo NULL -1
  """
}

process export_dada2tsv {

  publishDir "${params.outputdir}/dada2", mode: 'copy'
  publishDir "${params.outputdir}/allout", mode: 'copy'

  conda (params.enable_conda ? "${projectDir}/envs/pandas.yaml" : null)
  container "amancevice/pandas:1.3.5-slim"

  input:
    path dada_raw

  output:
    path 'asv.tab', emit: asvtab
    path 'repsep.fasta', emit: repsep

  """
    format_dada2.py ${dada_raw}
  """
}

//process clean_taxonomy_tsv {
//    publishDir "${params.outputdir}/allout", mode: 'copy'
//
//  container "amancevice/pandas:1.3.5-slim"
//
//  input:
//    path taxa_raw
//  
//  output:
//    path 'taxonomy.tsv', emit: taxtab
//}


/*
*  workflow
*/
workflow minflow {
    INPUT_CHECK ()
    ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
    ch_raw_long_reads  = INPUT_CHECK.out.raw_long_reads

    fastp(ch_raw_short_reads, params.trimfront1, params.trimfront2, params.trimtail1, params.trimtail2 )
    flash(fastp.out.reads)

    //  This does not work.
    // flash.out.reads.reduce([:])(map, tuple -> tuple[1]).subscribe onNext: { println it }, onComplete: { println 'Done' }
    chn_merge = flash.out.reads.flatMap(it -> it[1]).collect()
    dada2(chn_merge)
    export_dada2tsv(dada2.out.dada_raw_table)

    // Skip making tree for singularity since it produce error?
    //if (params.rootmafft == true) {
    // qiime2_roottree_mafft(export_dada2tsv.out.repsep)
    //}
    //else{
    //  qiime2_roottree(export_dada2tsv.out.repsep, export_dada2tsv.out.asvtab, params.roottree)
    //}

    classifier(export_dada2tsv.out.repsep, params.modeltax)
    picrust2(export_dada2tsv.out.repsep, export_dada2tsv.out.asvtab, "QIIME2", "test")
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
