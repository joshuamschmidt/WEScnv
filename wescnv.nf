#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.saveMode = 'copy'
params.inputFile = 'inputFile.txt'
params.batch = 'batch01'

// TODO: create input file channel
inputFile = file(params.inputFile)
lines = inputFile.readLines()

batch = params.batch
process cramCoverage {
    publishDir "$params.outdir/CoverageSummary", pattern: "*.summary.txt"

    label 'bamTasks'

    input:
    each line from lines

    output:
    file "${sample_id}.regions.bed.gz" into coverageOutChannel
    file "${sample_id}.mosdepth.summary.txt"

    script:
    list = line.split('\t')
    sample_id = list[0]
    input_cram = list[1]
    input_crai = list[2]
    """
    mosdepth --fasta \${ref_fasta} \
    --by \${target_bed} \
    --no-per-base \
    --mapq 25 \
    --threads $task.cpus \
    $sample_id \
    $input_cram
    """
}

process cramCounts {

    label 'bamTasks'

    input:
    each line from lines

    output:
    file "${sample_id}.cpt.bed.gz" into countsOutChannel

    script:
    list = line.split('\t')
    sample_id = list[0]
    input_cram = list[1]
    input_crai = list[2]
    """
    hts_nim_tools count-reads \
    --fasta \${ref_fasta} \
    --mapq 25 \
    --threads $task.cpus \
    \${target_bed} $input_cram \
    | sort -k1,1 -k2,2n \
    | gzip > "$sample_id".cpt.bed.gz
    """
}

countsOutChannel.into { aggregateCounts_ch; aggregateFpkm_ch }

process aggregateCoverage {
    publishDir "$params.outdir/CombinedCov/", pattern: "*coverage*"

    label 'combineTasks'

    input:
    file '*.regions.bed.gz' from coverageOutChannel.collect()

    output:
    "${batch}.coverage_MS-GC.GC5-DF-SD.bed.gz" into aggCountsOut_ch

    script:
    """
    countsToMatrix.py *.regions.bed.gz \
    --suffix .regions.bed.gz \
    --bed \${target_bed} \
    --merge-bed \
    | gzip > "$batch".coverage_MS-GC.GC5-DF-SD.bed.gz
    """
}


process aggregateCounts {
    publishDir "$params.outdir/CombinedCov/", pattern: "*counts*"

    label 'combineTasks'

    input:
    file '*.cpt.bed.gz' from aggregateCounts_ch.collect()

    output:
    "${batch}.counts_MS-GC.GC5-DF-SD.bed.gz"

    script:
    """
    countsToMatrix.py *.cpt.bed.gz \
    --suffix .cpt.bed.gz \
    --bed \${target_bed} \
    --merge-bed \
    | gzip > "$batch".counts_MS-GC.GC5-DF-SD.bed.gz
    """
}

process aggregateFpkm {
    publishDir "$params.outdir/CombinedCov", pattern: "*fkpm*"

    label 'combineTasks'

    input:
    file '*.cpt.bed.gz' from aggregateFpkm_ch.collect()

    output:
    "${batch}.fkpm_MS-GC.GC5-DF-SD.bed.gz"

    script:
    """
    countsToMatrix.py *.cpt.bed.gz \
    --suffix .cpt.bed.gz \
    --bed \${target_bed} \
    --fpkm \
    --merge-bed \
    | gzip > "$batch".fkpm_MS-GC.GC5-DF-SD.bed.gz
    """
}


/*
* Now, Let's take the counts file and determine which exons to remove due to consistently no data
*
*/


process filterCounts {

    input:
    file *counts*.gz from aggCountsOut_ch()

    output:
    "${batch}.counts_MS-GC.GC5-DF-SD.bed.gz"

}
