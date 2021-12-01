#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.saveMode = 'copy'
params.inputFile = 'inputFile.txt'

inputFile = file(params.inputFile)
lines = inputFile.readLines()

process foo {
    publishDir "$params.outdir/$sampleId", pattern: "*.txt"
    input:
    each line from lines

    output:
    file '*.txt'

    script:
    list = line.split('\t')
    sampleId = list[0]
    cram = list[1]
    crai = list[2]

    """
    echo $sampleId $cram $crai > test.txt
    """
}
