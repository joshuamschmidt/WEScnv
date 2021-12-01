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
    publishDir "$params.outdir/CoverageSummary", mode: "${saveMode}", pattern: "*.summary.txt"
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
    label 'combineTasks'

    publishDir "$params.outdir/CombinedCov", mode: "${saveMode}", pattern: "*coverage_MS-GC.GC5-DF-SD.bed.gz"

    input:
    file '*.regions.bed.gz' from coverageOutChannel.collect()

    output:
    "$batch.coverage_MS-GC.GC5-DF-SD.bed.gz"

    script:
    """
    set -euo pipefail
    countsToMatrix.py *.regions.bed.gz \
    --suffix .regions.bed.gz \
    --bed \${target_bed} \
    --merge-bed \
    | gzip > "$batch".coverage_MS-GC.GC5-DF-SD.bed.gz
    """
}


process aggregateCounts {
    label 'combineTasks'

    publishDir "$params.outdir/CombinedCov", mode: "${saveMode}", pattern: "*counts_MS-GC.GC5-DF-SD.bed.gz"
    input:
    file '*.cpt.bed.gz' from aggregateCounts_ch.collect()

    output:
    "$batch.counts_MS-GC.GC5-DF-SD.bed.gz"

    script:
    """
    set -euo pipefail
    countsToMatrix.py *.cpt.bed.gz \
    --suffix .cpt.bed.gz \
    --bed \${target_bed} \
    --merge-bed \
    | gzip > "$batch".counts_MS-GC.GC5-DF-SD.bed.gz
    """
}

process aggregateFpkm {
    label 'combineTasks'

    publishDir "$params.outdir/CombinedCov", mode: "${saveMode}", pattern: "*fpkm_MS-GC.GC5-DF-SD.bed.gz"

    input:
    file '*.cpt.bed.gz' from aggregateFpkm_ch.collect()

    output:
    "$batch.fkpm_MS-GC.GC5-DF-SD.bed.gz"

    script:
    """
    set -euo pipefail
    countsToMatrix.py *.cpt.bed.gz \
    --suffix .cpt.bed.gz \
    --bed \${target_bed} \
    --fpkm \
    --merge-bed \
    | gzip > "$batch".fkpm_MS-GC.GC5-DF-SD.bed.gz
    """
}
