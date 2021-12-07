#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.saveMode = 'copy'
params.inputFile = 'inputFile.txt'
params.batch = 'batch01'
batch=params.batch

Channel
    .fromPath(params.inputFile)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.sample_id, file(row.input_cram), file(row.input_crai)) }
    .set { samples_ch }

samples_ch.into { coverageInChannel; countsInChannel }

process cramCoverage {
    publishDir "$params.outdir/CoverageSummary", pattern: "*.summary.txt"

    label 'bamTasks'

    input:
    set sample_id, file(input_cram), file(input_crai) from coverageInChannel

    output:
    file "*regions.bed.gz" into coverageOutChannel
    file "${sample_id}.mosdepth.summary.txt" into coverageSummaryChannel

    script:
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
    set sample_id, file(input_cram), file(input_crai) from countsInChannel

    output:
    file "*.cpt.bed.gz" into countsOutChannel

    script:
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
    file input_files from coverageOutChannel.collect()
    //file "*" from coverageOutChannel.collect()

    output:
    file "${batch}.coverage_MS-GC.GC5-DF-SD.bed.gz"

    script:
    """
    countsToMatrix.py $input_files \
    --suffix ".regions.bed.gz" \
    --bed \${target_bed} \
    --merge-bed \
    | gzip > "$batch".coverage_MS-GC.GC5-DF-SD.bed.gz
    """
}


process aggregateCounts {
    publishDir "$params.outdir/CombinedCov/", pattern: "*counts*"

    label 'combineTasks'

    input:
    file input_files from aggregateCounts_ch.collect()

    output:
    file "${batch}.counts_MS-GC.GC5-DF-SD.bed.gz"

    script:
    """
    countsToMatrix.py $input_files \
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
    file input_files from aggregateFpkm_ch.collect()

    output:
    file "${batch}.fkpm_MS-GC.GC5-DF-SD.bed.gz"

    script:
    """
    countsToMatrix.py $input_files \
    --suffix .cpt.bed.gz \
    --bed \${target_bed} \
    --fpkm \
    --merge-bed \
    | gzip > "$batch".fkpm_MS-GC.GC5-DF-SD.bed.gz
    """
}


process assignBioSex {
    publishDir "$params.outdir/CombinedCov/", pattern: "*Assignment.txt"

    label 'combineTasks'

    input:
    file input_files from coverageSummaryChannel.collect()

    output:
    file "${batch}_bioSex-Assignment.txt"

    script:
    """
    assign_bio_sex.py $input_files \
    --suffix ".regions.bed.gz" > "$batch"_bioSex-Assignment.txt
    """
}



/*
* Now, Let's take the fkpm data to find clusters of samples (define) sub-batches
* and define which exons should be filtered from the analysis.
* outputs plots of clusterings...
* compare to ExomeDepth method?
*/


/*process filterCounts {

    input:
    file '*counts*.gz' from aggCountsOut_ch()

    output:
    "${batch}.counts_MS-GC.GC5-DF-SD.bed.gz"

}
*/
