#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.saveMode = 'copy'
params.inputFile = 'inputFile.txt'


Channel
    .fromPath(params.inputFile)
    .splitCsv(header:true)
    .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    .set { samples_ch }

process foo {
    input:
    set sampleId, file(read1), file(read2) from samples_ch

    script:
    """
    echo "your_command --sample $sampleId --reads $read1 $read2" > $sampleId.txt
    """
}
