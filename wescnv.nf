#!/usr/bin/env nextflow

// default params
params.outputDir = 'run/'
params.saveMode = 'copy'
params.inputFile = 'inputFile.txt'

inputFile = file(params.inputFile)
lines = inputFile.readLines()

process cramCoverage {
    publishDir "$params.outdir/CoverageSummary", pattern: "*.summary.txt"
    label 'bamTasks'

    input:
    each line from lines

    output:
    file "${sample_id}.regions.bed.gz" into coverageOutChannel

    script:
    list = line.split('\t')
    sample_id = list[0]
    input_cram = list[1]
    input_crai = list[2]
    output_coverage_filename = "${sample_id}.regions.bed.gz"
    coverage_summary_filename = "${sample_id}.mosdepth.summary.txt"
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
    output_counts_filename = "${sample_id}.cpt.bed.gz"
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

coverageOutChannel
        .view()


countsOutChannel
        .view()

/*
counts = Channel.fromPath( '/some/path/*.fa' )

process blastThemAll {
  input:
  file counts

  "blastp -query $counts -db nr"

}

*/

/*
process aggregateCounts {

    input:
    file '*.CPT.txt.gz' from nimtoolsCountsOutChannel.collect()

    output:
    set val(ius), file("${recalibration_report_filename}") into BaseRecalOut

    script:
    sequence_group_interval = sequence_group_interval_all.get(i).toString()
    if ( i < 10 ) {
        recalibration_report_filename = "${ius}.0${i}.recal_data.csv"
    } else {
        recalibration_report_filename = "${ius}.${i}.recal_data.csv"
    }
    """
    gatk --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -XX:+PrintFlagsFinal \
        -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
        -Xloggc:gc_log.log -Xms4000m" \
        BaseRecalibrator \
        -R \${ref_fasta} \
        -I ${input_bam} \
        --use-original-qualities \
        -O ${recalibration_report_filename} \
        -known-sites \${dbSNP_vcf} \
        -known-sites \${known_indels_sites_Gold} \
        -known-sites \${known_indels_sites_assembly38} \
        -L \$(echo ${sequence_group_interval} | sed 's/ / -L /g')
    """
}
*/

