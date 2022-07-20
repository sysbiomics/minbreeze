/*
 *  DAF pipeline
 *
 */

nextflow.enable.dsl=2

include { minflow } from './workflows/minflow.nf'

workflow NFCORE_MINFLOW {
    minflow ()
}

workflow {
    NFCORE_MINFLOW ()
}
