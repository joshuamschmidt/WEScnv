#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.saveMode = 'copy'
params.inputFile = 'inputFile.txt'


Channel
    .fromPath(params.inputFile)
    .splitCsv(header:true, sep:"\t")
    .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
    .set { samples_ch }

process foo {
    publishDir "$params.outdir/$sampleId"
    input:
    set sampleId, file(read1), file(read2) from samples_ch

    output:
    file '*.txt'

    script:
    """
    echo $sampleId $read1 $read2 > test.txt
    """
}