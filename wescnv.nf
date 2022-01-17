#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.inputFile = 'inputFile.txt'
params.batch = 'batch01'
params.reference_fasta = 'ref.fasta'
params.reference_fasta_index = 'ref.fai'
params.target_bed = 'targets.bed'
params.target_cov_txt = 'targets.txt'
params.target_picard = ''
params.bait_bed = 'baits.bed'
params.bait_picard = ''
batch=params.batch
reference_fasta=params.reference_fasta
reference_fasta_index=params.reference_fasta_index
target_bed=params.target_bed
target_cov_txt=params.target_cov_txt
target_picard_list=params.target_picard
bait_bed=params.bait_bed
bait_picard_list=params.bait_picard


Channel
    .fromPath(params.inputFile)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.sample_id, file(row.input_cram), file(row.input_crai)) }
    .set { samples_ch }

samples_ch.into { coverageInChannel; countsInChannel; hsMetricsInChannel; isMetricsInChannel}

process cramCoverage {
    publishDir "$params.outdir/CoverageSummary", pattern: "*.summary.txt"

    label 'bamTasks'

    input:
    tuple val(sample_id), path(input_cram), path(input_crai) from coverageInChannel

    output:
    path "*regions.bed.gz" into coverageOutChannel
    path "${sample_id}.mosdepth.summary.txt" into coverageSummaryChannel

    script:
    """
    cp $reference_fasta_index .
    mosdepth --fasta $reference_fasta \
    --by $target_bed \
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
    tuple val(sample_id), path(input_cram), path(input_crai) from countsInChannel

    output:
    path "*.cpt.bed.gz" into countsOutChannel

    script:
    """
    cp $reference_fasta_index .
    hts_nim_tools count-reads \
    --fasta $reference_fasta \
    --mapq 25 \
    --threads $task.cpus \
    $target_bed $input_cram \
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

    output:
    file "${batch}.coverage.bed.gz"

    script:
    """
    countsToMatrix.py $input_files \
    --suffix ".regions.bed.gz" \
    --bed $target_cov_txt \
    --merge-bed \
    | gzip > "$batch".coverage.bed.gz
    """
}

process aggregateCounts {
    publishDir "$params.outdir/CombinedCov/", pattern: "*counts*"

    label 'combineTasks'

    input:
    file input_files from aggregateCounts_ch.collect()

    output:
    file "${batch}.counts.bed.gz"

    script:
    """
    countsToMatrix.py $input_files \
    --suffix .cpt.bed.gz \
    --bed $target_cov_txt \
    --merge-bed \
    | gzip > "$batch".counts.bed.gz
    """
}

process aggregateFpkm {
    publishDir "$params.outdir/CombinedCov", pattern: "*fkpm*"

    label 'combineTasks'

    input:
    file input_files from aggregateFpkm_ch.collect()

    output:
    file "${batch}.fkpm.bed.gz"

    script:
    """
    countsToMatrix.py $input_files \
    --suffix .cpt.bed.gz \
    --bed $target_cov_txt \
    --fpkm \
    --merge-bed \
    | gzip > "$batch".fkpm.bed.gz
    """
}

process assignBioSex {
    publishDir "$params.outdir/AssignedSex/", pattern: "*Assignment.txt"

    label 'combineTasks'

    input:
    file input_files from coverageSummaryChannel.collect()

    output:
    file "${batch}_bioSex-Assignment.txt"

    script:
    """
    assign_bio_sex.py $input_files \
    --suffix ".mosdepth.summary.txt" > "$batch"_bioSex-Assignment.txt
    """
}


process collectHSMetrics {
    publishDir "$params.outdir/HSmetrics/", pattern: "*hs_metrics.txt"
    label 'picardMetrics'

    input:

    tuple val(sample_id), path(input_cram), path(input_crai) from hsMetricsInChannel

    output:

    tuple val(sample_id), path("*hs_metrics.txt") into HSMetricsOuts

    script:
    """
    cp $reference_fasta_index .
    gatk --java-options "-Xmx8g" CollectHsMetrics \
      I=$input_cram \
      O="$sample_id"_hs_metrics.txt \
      R=$reference_fasta \
      BAIT_INTERVALS=$bait_picard_list \
      TARGET_INTERVALS=$target_picard_list
    """
}

process collectISMetrics {
    publishDir "$params.outdir/ISmetrics/", pattern: "*is_metrics*"
    label 'picardMetrics'

    input:

    tuple val(sample_id), path(input_cram), path(input_crai) from isMetricsInChannel

    output:

    tuple val(sample_id), path("*is_metrics.txt") into ISMetricsOuts

    script:
    """
    cp $reference_fasta_index .
    gatk --java-options "-Xmx8g" CollectInsertSizeMetrics \
      R=$reference_fasta \
      I=$input_cram \
      O="$sample_id"_is_metrics.txt \
      H="$sample_id"_is_metrics.pdf
    """
}

mergedMetricsInChannel = HSMetricsOuts.join(ISMetricsOuts, failOnDuplicate: true, failOnMismatch: true)

process MergeMetrics{
    publishDir "$params.outdir/Mergedmetrics/", pattern: "*_mergedMetrics.txt"
    label 'script'
    input:

    tuple val(sample_id), path(hs_metrics), path(is_metrics) from mergedMetricsInChannel

    output:

    path("*_mergedMetrics.txt") into defineClustersInChannel

    script:
    """
    #!/usr/bin/env bash
    set -eo pipefail
    echo $sample_id > tmp_name
    cut -f9,32,46,48,52,63,64 $hs_metrics | sed '8q;d' > tmp_hs
    cut -f1,3,6,7 $is_metrics | sed '8q;d'> tmp_is
    paste tmp_name tmp_hs tmp_is > "$sample_id"_mergedMetrics.txt
    """
}

/*
process defineProcessGroups {

    input:
    file input_files from defineClustersInChannel.collect()

}

*/


/*process filterCounts {

    input:
    file '*counts*.gz' from aggCountsOut_ch()

    output:
    "${batch}.counts_MS-GC.GC5-DF-SD.bed.gz"

}
*/
