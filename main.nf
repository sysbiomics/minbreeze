/*
*  DAF pipeline
*
*/

// For loading config from yaml
import java.nio.file.Files
import java.nio.file.Paths
import org.yaml.snakeyaml.Yaml


DAF_VERSION="0.00a"

def helpMessage() {
  log.info"""
  --inputdir     Path to input data
  --outputdir    Path to output data
  --manifest     Path to manifest
  --config       Path to config
  """.stripIndent()
}

/* Configuration
*/

println("${params.getClass()}")
// Show help message
if (params.help){
    helpMessage()
    exit 0
}

// Read config files
pmatcher = "*_R{1,2}.fastq.gz"
matcher = params.inputdir + "/" + pmatcher

// Start moving config to config file
conf = [:]

// Manifest (check if == NONE)
def extconf = (Map)new Yaml().load(Files.newInputStream(Paths.get("./param.conf")))

// Primers
conf.fwdprimer = "CCTACGGGNGGCWGCAG"
conf.revprimer = "GACTACHVGGGTATCTAATCC"
conf.fwdprimerlen = extconf["fwdprimerlen"]
conf.revprimerlen = extconf["revprimerlen"]
conf.model = extconf["model"]
// Denoiser
// Classifier
conf.modelabs = new File(conf.model).getCanonicalPath();
// Root tree for fragment insertion
conf.roottree="./res/sepp-refs-gg-13-8.qza"
conf.roottreeabs=new File(conf.roottree).getCanonicalPath();




/* End of configuration
*/


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

  publishDir "${params.outputdir}", mode: 'copy'
  input:
  set val(pair_id), file(reads) from ch_read_pairs

  output:
  file("trimmed/${reads[0]}") into ch_dada2forReads
  file("trimmed/${reads[1]}") into ch_dada2revReads
  file("trimmed/*.log")

  """
  mkdir -p t_trimleft
  mkdir -p trimmed
  echo ${pair_id}
  seqtk trimfq -b ${conf.fwdprimerlen} ${reads[0]} | gzip > t_trimleft/${reads[0]}
  seqtk trimfq -b ${conf.revprimerlen} ${reads[1]} | gzip > t_trimleft/${reads[1]}

  bbduk.sh -Xmx1g in=t_trimleft/${reads[0]} in2=t_trimleft/${reads[1]} \
  out1=trimmed/${reads[0]} out2=trimmed/${reads[1]} \
  qtrim=r trimq=15 \
  minlength=150 stats=trimmed/stat_${pair_id}.log \
  2> trimmed/run_${pair_id}.log
  echo "Deterministic trim with ${conf.fwdprimerlen} ${conf.revprimerlen}" > trimmed/seqtk.log
  """
}

//Trim primer/adapter with matching
// Current set to false
// Maybe I should use this AFTER det trim to make sure that everything is clean

if (false){
    literals=conf.fwdprimer + "," + conf.revprimer
    process heutrim {

    publishDir "${params.outputdir}", mode: 'copy'
    input:
    set val(pair_id), file(reads) from ch_read_pairs

    output:
    file("trimmed/${reads[0]}") into ch_dada2forReads
    file("trimmed/${reads[1]}") into ch_dada2revReads
    file("trimmed/*.log")

    """
    mkdir -p trimmed
    echo ${pair_id}
    bbduk.sh -Xmx1g in=${reads[0]} in2=${reads[1]} out1=trimmed/${reads[0]} \
      out2=trimmed/${reads[1]} \
      literal=${literals} \
      ktrim=l k=13 mink=6 hdist=1 qtrim=r trimq=15 \
      minlength=150 restrictleft=25 copyundefined=T stats=trimmed/stat_${pair_id}.log \
      2> trimmed/run_${pair_id}.log
    """
    }
}


/* Denoise
*/
process dada2 {
  label 'big_cpu'

  publishDir "${params.outputdir}/dada2", mode: "copy"
  input:
    file freads from ch_dada2forReads.collect()
    file rreads from ch_dada2revReads.collect()

  output:
    file("dadaraw.tsv") into ch_dada2ext
    file("track.tsv")
    file("qualplotF.pdf")
    file("qualplotR.pdf")
  """
  mkdir -p fwd
  mkdir -p rev
  mv ${freads} fwd
  mv ${rreads} rev
  run_dada_paired.R fwd rev dadaraw.tsv track.tsv ff rf 0 0 0 0 2.0 2.0 2 consensus 1.0 ${task.cpus} 1000000
  """
}

process dada2ext {
  publishDir "${params.outputdir}/dada2", mode: "copy"

  input:
    file rawasv from ch_dada2ext

  output:
    file("asv.tab")
    file("repsep.fasta") into ch_qiime2_roottree
    file("repsep.fasta") into ch_qiime2_bayes

  """
  format_dada2.py ${rawasv} ${conf.manifest} # Produce asv.tab and repsep.fasta
  """
}

/*
  Building phylogenetic tree
*/

process qiime2_roottree {

  label 'mod_cpu'
  publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'
  input:
    file repsep_fasta from ch_qiime2_roottree
  output:
    file("rooted-tree.qza")
    file("rooted-tree.nwk")

  """
  qiime tools import --input-path ${repsep_fasta} \
    --output-path sequences.qza \
    --type 'FeatureData[Sequence]'

  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences sequences.qza \
    --o-alignment aligned-rep-seqs.qza \
    --o-masked-alignment masked-aligned-rep-seqs.qza \
    --o-tree unrooted-tree.qza \
    --o-rooted-tree rooted-tree.qza

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

  //#qiime fragment-insertion sepp \
  //#  --i-representative-sequences ./sequences.qza \
  //#  --i-reference-database ${conf.roottreeabs} \
  //#  --o-tree ./rooted-tree.qza \
  //#  --o-placements ./tree_placements.qza \
  //#  --p-threads 1  # update to a higher number if you can

if (false) {
process qiime2_blast {

  label 'mod_cpu'
  publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'


  input:
    file repsep_fasta from ch_qiime2_bayes

  output:
    file("taxonomy.qza")
    file("taxonomy.tsv")

  """
  model=${conf.modelabs}

  qiime tools import --input-path ${repsep_fasta} \
    --output-path sequences.qza \
    --type 'FeatureData[Sequence]'

  qiime feature-classifier classify-consensus-blast \
    --i-query sequences.qza \
    --i-reference-reads ${baseDir}/res/gg99_seq.qza \
    --i-reference-taxonomy ${baseDir}/res/gg99_tax.qza \
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
}
process qiime2_bayes {

  label 'mod_cpu'
  publishDir "${params.outputdir}/qiime2_analysis", mode: 'copy'

  input:
    file repsep_fasta from ch_qiime2_bayes

  output:
    file("taxonomy.qza")
    file("taxonomy.tsv")

  """
  model=${conf.modelabs}

  qiime tools import --input-path ${repsep_fasta} \
    --output-path sequences.qza \
    --type 'FeatureData[Sequence]'

  qiime feature-classifier classify-sklearn \
    --i-classifier \${model} \
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

process param_report {

  publishDir "${params.outputdir}"

  shell:
  '''
  echo "Export all param used into files"
  '''
}

process report_gen {
  publishDir "${params.outputdir}"

  shell:
  '''
  echo "Export all param used into files"
  '''

}
