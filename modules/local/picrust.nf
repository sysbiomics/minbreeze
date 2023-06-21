//MIT License
//
//Copyright (c) Daniel Straub, Alexander Peltzer
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

process picrust2 {
    tag "${seq},${abund}"
    label 'process_medium'

    conda "bioconda::picrust2=2.4.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picrust2:2.4.2--pyhdfd78af_0' :
        'quay.io/biocontainers/picrust2:2.4.2--pyhdfd78af_0' }"

    publishDir "${params.outputdir}/picrust2_analysis", mode: 'copy'

    input:
      path(seq)
      path(abund)
      val(source)
      val(message)

    output:
      path("picrust2_analysis/*") , emit: outfolder
      path("*_descrip.tsv"), emit: pathways
      path "versions.yml"  , emit: versions
      path "*.args.txt"    , emit: args
      path "${message}.txt"

    script:
    def args = task.ext.args ?: ''
    """
    #If input is QIIME2 file, than (1) the first line and (2) the first character (#) of the second line need to be removed
    if [ "$source" == 'QIIME2' ]
    then
        tail -n +2 "$abund" > "${abund}.tmp" && mv "${abund}.tmp" "$abund"
    fi

    picrust2_pipeline.py \\
        $args \\
        -s $seq \\
        -i $abund \\
        -o picrust2_analysis \\
        -p $task.cpus \\
        --in_traits EC,KO \\
        --verbose

    #Add descriptions to identifiers
    add_descriptions.py -i picrust2_analysis/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \
                    -o EC_pred_metagenome_unstrat_descrip.tsv
    add_descriptions.py -i picrust2_analysis/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \
                    -o KO_pred_metagenome_unstrat_descrip.tsv
    add_descriptions.py -i picrust2_analysis/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \
                    -o METACYC_path_abun_unstrat_descrip.tsv

    echo "$message" > "${message}.txt"
    echo -e "picrust\t$args" > "picrust.args.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        picrust2: \$( picrust2_pipeline.py -v | sed -e "s/picrust2_pipeline.py //g" )
    END_VERSIONS
    """
}
