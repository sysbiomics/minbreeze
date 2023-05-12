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
params.fprimer = 'CCTAYGGGRBGCASCAG'
params.rprimer = 'GGACTACNNGGGTATCTAAT'
// Or
params.fcut = 10
params.rcut = 15
/* Skip some analysis
*/
params.SKIP_NWKTREE = false
params.SKIP_SEQCLAS = false
params.SKIP_FUNCLAS = false

params.db = 'silva' // or gg
params.flash_merge = true
params.sepp_tree = false

// For BLAST
//  blastrefseq = "s3://qiime2-data/2021.8/common/"
//  blastreftax = "s3://qiime2-data/2021.8/common/"

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
// Module
//

include { picrust2                                   } from '../modules/local/picrust'
include { dada2_single; dada2_pair; export_dada2tsv  } from '../modules/local/dada2'

//
// Workflow
//

include { INPUT_CHECK                                 } from '../subworkflows/local/input_check'
include { qiime2_bayes                                } from '../subworkflows/local/classify'
include { classify_reads                              } from '../subworkflows/local/classify'
include { qiime2_roottree_mafft; qiime2_roottree_sepp } from '../subworkflows/local/phylotree'


/*
 * Trimming primer and quality. We don't use DADA2 for this because user 
 * might want to try other trimming solution instead.
*/
process fastp {

  conda (params.enable_conda ? "${projectDir}/envs/minflow.yaml" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'quay.io/biocontainers/fastp:0.23.2--h5f740d0_3' :
      'quay.io/biocontainers/fastp:0.23.2--h5f740d0_3' }"

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

// Primer removal*
/*process cutadapt {
}
*/

/* Merge read
TODO: Make output compatible with later
*/
process flash {

  conda (params.enable_conda ? "${projectDir}/envs/minflow.yaml" : null)
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'quay.io/biocontainers/flash:1.2.11--hed695b0_5' :
      'quay.io/biocontainers/flash:1.2.11--hed695b0_5' }"

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


/*
*  workflow
*/
workflow minflow {
    INPUT_CHECK ()
    ch_raw_short_reads = INPUT_CHECK.out.raw_short_reads
    // ch_raw_long_reads  = INPUT_CHECK.out.raw_long_reads

    fastp(ch_raw_short_reads, params.trimfront1, params.trimfront2, params.trimtail1, params.trimtail2 )

    if (params.flash_merge) { // Merge and use dada2 single
      flash(fastp.out.reads)
      chn_merge = flash.out.reads.flatMap(it -> it[1]).collect()
      dada2_single(chn_merge)
      export_dada2tsv(dada2_single.out.dada_raw_table)
    }
    else {  // Let's DADA2 handle everythings
      fastp.out.reads.multiMap { it ->
            fwd: it[1][0]
            rev: it[1][1]}.set { rawseqs }
      fwdall = rawseqs.fwd.collect()
      revall = rawseqs.rev.collect()

      dada2_pair(fwdall, revall)
      export_dada2tsv(dada2_pair.out.dada_raw_table)
    }

    if (! params.SKIP_NWKTREE) {
        if (params.sepp_tree) {
            qiime2_roottree_sepp(export_dada2tsv.out.repsep, export_dada2tsv.out.asvtab, params.sepp_tree)
        }
        else {
            qiime2_roottree_mafft(export_dada2tsv.out.repsep)
        }
    }

    if (! params.SKIP_SEQCLAS) {
        classify_reads(export_dada2tsv.out.repsep, params.modeltax, null, null)
        // qiime2_bayes(export_dada2tsv.out.repsep, params.modeltax)
    }

    if (! params.SKIP_FUNCLAS) {
        picrust2(export_dada2tsv.out.repsep, export_dada2tsv.out.asvtab, export_dada2tsv.out.source, "qiime2")
    }
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

workflow.onError {
    Send.sendMessage('OOP', 'localhost', 12020)
    // println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}

// vi: ft=groovy
