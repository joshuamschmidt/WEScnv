
// aln is either a cram or bam file, and its proper index


Channel
    .fromPath('/Users/joshuaschmidt/testlist.txt')
    .splitCsv(header:true, sep:'\t')
    .map { create_aln_channel(it) }
    .view()


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
    
