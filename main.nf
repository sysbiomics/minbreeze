/*
*  DAF pipeline
*
*/

// For loading config from yaml
import java.nio.file.Files
import java.nio.file.Paths

nextflow.preview.dsl=2

def helpMessage() {
  log.info"""
  --inputdir     Path to input data
  --outputdir    Path to output data
  --config       Path to config
  """.stripIndent()
}

// Help message
if (params.help){
    helpMessage()
    exit 0
}
/* Configuration
*/
// Read config files
pmatcher = params.matcher
matcher = params.inputdir + "/" + pmatcher

// Cleaning config in case something does not get used,
// For example, QIIME2_blast does not use bayes model.
def validateParams(){
}

params.modelabs = new File(params.model).getCanonicalPath();
// Root tree for fragment insertion
params.roottreeabs= new File(params.modeltree).getCanonicalPath();
// params.manifestabs= new File(params.manifest).getCanonicalPath();
// BLAST or vsearch
params.refseqabs = new File(params.refseq).getCanonicalPath();
params.reftaxabs = new File(params.reftax).getCanonicalPath();

/* End of configuration
*/

process fastqc {

  publishDir "${params.outputdir}/fastqc", mode: 'copy'

  input:
    tuple val(pair_id), file(reads)
  output:
    file "*.html"
    file "*.zip"

  shell:
  """
  fastqc ${reads[0]}
  fastqc ${reads[1]}
  """
}

/*
  Trimming primer and quality. We don't use DADA2 for this because we *might*
  want to try other pipeline instead.
*/

/*  Determistic trimming, trim the front with set length.
*/
process dettrim {

  publishDir "${params.outputdir}/trimmed", mode: 'link'

  input:
    tuple val(pair_id), file(reads)
    val(truncleft)
    val(truncright)
  output:
    path "trimmed/${reads[0]}", emit: fwd_reads
    path "trimmed/${reads[1]}", emit: rev_reads
    path "trimmed/*.log"

  """
  mkdir -p t_trimleft
  mkdir -p trimmed
  echo ${pair_id}
  seqtk trimfq -b ${truncleft} ${reads[0]} | gzip > t_trimleft/${reads[0]}
  seqtk trimfq -b ${truncright} ${reads[1]} | gzip > t_trimleft/${reads[1]}

  bbduk.sh -Xmx1g in=t_trimleft/${reads[0]} in2=t_trimleft/${reads[1]} \
  out1=trimmed/${reads[0]} out2=trimmed/${reads[1]} \
  qtrim=r trimq=15 \
  minlength=150 stats=trimmed/stat_${pair_id}.log \
  2> trimmed/run_${pair_id}.log
  echo "Deterministic trim with ${params.fwdprimerlen} ${params.revprimerlen}" > trimmed/seqtk.log
  """
}

/* Denoise
*/
process dada2 {
  label 'big_cpu', 'dada2'

  publishDir "${params.outputdir}/dada2", mode: "copy"
  input:
    path freads
    path rreads
    val fwdlen
    val revlen
    val trimleft
    val trimright

  output:
    path "dadaraw.tsv", emit: dada_raw_table
    path "track.tsv", emit: track
    path "qualplotF.pdf"
    path "qualplotR.pdf"
  """
  mkdir -p fwd
  mkdir -p rev
  mv ${freads} fwd
  mv ${rreads} rev
  run_dada_paired.R fwd rev dadaraw.tsv track.tsv ff rf ${fwdlen} ${revlen} ${trimleft} ${trimright} 2.0 2.0 2 ${params.dada.chimera_alg} ${params.dada.chimera_fol} ${task.cpus} 1000000 ${params.dada.pool}
  """
}

process dada2ext {

  publishDir "${params.outputdir}/dada2", mode: "copy"

  input:
    path rawasv

  output:
    path "asv.tab"
    path "repsep.fasta", emit: repsep

  script:
  if(params.manifestabs)
    """
    format_dada2.py ${rawasv} ${params.manifestabs}
    """
  else
    """
    format_dada2.py ${rawasv}
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
  output:
    path "rooted-tree.qza"
    path "rooted-tree.nwk"

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
    --i-reference-database ${params.roottreeabs} \
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

include './libs/classify.nf'

/*
*  Final workflow
*/

// Single (forward)

/*
include './libs/single.nf'
workflow {
chan = Channel
    .fromPath(matcher)
    .ifEmpty { exit 1, "Cannot find anything here chief"}
fastqc_single(chan)
dettrim_single(chan)
publish:
    fastqc.out to: "${params.outputdir}/fastqc", mode: 'copy'
    dettrim.out to: "${params.outputdir}/trimmed", mode: 'link'
}
*/

/*
  if (params.qtrim==true){
    dettrim(chan, params.fwdprimerlen, params.revprimerlen)
    dada2(dettrim.out.fwd_reads.collect(), dettrim.out.rev_reads.collect(), params.fwdlen, params.revlen, params.fwdprimerlen, params.revprimerlen)
  } else {
    dada2(something.fwd.collect(), something.rev.collect(), params.fwdlen, params.revlen, params.fwdprimerlen, params.revprimerlen)
  }
  */

workflow {
  chan = Channel
    .fromFilePairs(matcher)
    .ifEmpty { exit 1, "Cannot find anything here chief"}
  fastqc(chan)
  if(params.qtrim==true){
    dettrim(chan, params.fwdprimerlen, params.revprimerlen)
    dada2(dettrim.out.fwd_reads.collect(), dettrim.out.rev_reads.collect(), params.fwdlen, params.revlen, params.fwdprimerlen, params.revprimerlen)
  } else {
        chan.multiMap { it ->
            fwd: it[1][0]
            rev: it[1][1]
        }.set { something }
        dada2(something.fwd.collect(), something.rev.collect(), params.fwdlen, params.revlen, params.fwdprimerlen, params.revprimerlen)
  }
  
  dada2ext(dada2.out.dada_raw_table)
  qiime2_roottree(dada2ext.out.repsep)
  qiime2_bayes(dada2ext.out.repsep, params.model)
  //qiime2_blast(dada2ext.out.repsep, params.refseqabs, params.reftaxabs)
}


// Export options use,
// For now, just export everything
import groovy.json.JsonBuilder
workflow.onComplete {
    // serialize
    def builder = new JsonBuilder()
    builder(params)
    def output_f = new File("${params.outputdir}/pipeline.txt")
    output_f.withWriter {w -> w << builder.toString()}
}
