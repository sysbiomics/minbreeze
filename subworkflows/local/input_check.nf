//
// Check input samplesheet and get read channels
//


params.options = [:]

def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

workflow INPUT_CHECK {
    main:
    if(hasExtension(params.input, "csv")){
        // extracts read files from samplesheet CSV and distribute into channels
        Channel
            .from(file(params.input))
            .splitCsv(header: true)
            .map { row ->
                    if (row.size() == 4) {
                        def id = row.sample
                        def group = row.group
                        def sr1 = row.fastq_1 ? file(row.fastq_1, checkIfExists: true) : false
                        def sr2 = row.fastq_2 ? file(row.fastq_2, checkIfExists: true) : false
                        // Check if given combination is valid
                        if (!sr1) exit 1, "Invalid input samplesheet: fastq_1 can not be empty."
                        if (!sr2 && !params.single_end) exit 1, "Invalid input samplesheet: single-end short reads provided, but command line parameter `--single_end` is false. Note that either only single-end or only paired-end reads must provided."
                        if (sr2 && params.single_end) exit 1, "Invalid input samplesheet: paired-end short reads provided, but command line parameter `--single_end` is true. Note that either only single-end or only paired-end reads must provided."
                        return [ id, group, sr1, sr2 ]
                    } else {
                        exit 1, "Input samplesheet contains row with ${row.size()} column(s). Expects 4."
                    }
                }
            .set { ch_input_rows }
        // separate short and long reads
        ch_input_rows
            .map { id, group, sr1, sr2 ->
                        def meta = [:]
                        meta.id           = id
                        meta.group        = group
                        meta.single_end   = params.single_end
                        if (params.single_end)
                            return [ meta, [ sr1] ]
                        else
                            return [ meta, [ sr1, sr2 ] ]
                }
            .set { ch_raw_short_reads }
    } else {
        Channel
            .fromFilePairs(params.input, size: params.single_end ? 1 : 2)
            .ifEmpty { exit 1, "Cannot find any reads matching: ${params.input}\nNB: Path needs to be enclosed in quotes!\nIf this is single-end data, please specify --single_end on the command line." }
            .map { row ->
                        def meta = [:]
                        meta.id           = row[0]
                        meta.group        = 0
                        meta.single_end   = params.single_end
                        return [ meta, row[1] ]
                }
            .set { ch_raw_short_reads }
        ch_input_rows = Channel.empty()
        ch_raw_long_reads = Channel.empty()
    }

    // Ensure sample IDs are unique
    ch_input_rows
        .map { id, group, sr1, sr2 -> id }
        .toList()
        .map { ids -> if( ids.size() != ids.unique().size() ) {exit 1, "ERROR: input samplesheet contains duplicated sample IDs!" } }

    emit:
    raw_short_reads = ch_raw_short_reads
}
