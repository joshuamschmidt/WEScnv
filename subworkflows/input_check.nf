//
// Check input samplesheet and get gsfile channels
// TODO
// Check format of samplesheet. At the moment it just queues gsfile channel
//
// Check input samplesheet and get aln/index channel
//


include { SAMPLESHEET_CHECK } from '../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:'\t' )
        .map { create_aln_channel(it) }
        .set { alignments }


    emit:
    alignments // channel: [ val(meta), [ aln, index ] ]
    // TODO versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ aln, index ] ]
// following example from https://github.com/nf-core/rnaseq/blob/master/subworkflows/local/input_check.nf

def create_aln_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id  = row.sample_id

    // add paths of the aln file and its index
    def aln_meta = []
    if (!file(row.sample_aln).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Alignment file (bam/cram) file does not exist!\n${row.sample_aln}"
    }
    if (file(row.sample_index).exists()) {
        aln_meta = [ meta, [ file(row.sample_aln), file(row.sample_index) ] ]
    } else {
        exit 1, "ERROR: Please check input samplesheet -> Index file does not exist!\n${row.sample_index}"
        }
    return aln_meta
   }


