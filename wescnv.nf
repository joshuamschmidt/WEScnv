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
params.cnvKitTarget = ''
params.cnvKitAntiTarget = ''
batch=params.batch
reference_fasta=params.reference_fasta
reference_fasta_index=params.reference_fasta_index
target_bed=params.target_bed
target_cov_txt=params.target_cov_txt
target_picard_list=params.target_picard
bait_bed=params.bait_bed
bait_picard_list=params.bait_picard
cnvkit_target_bed=params.cnvKitTarget
cnvkit_antitarget_bed=params.cnvKitAntiTarget

Channel
    .fromPath(params.inputFile)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.sample_id, file(row.input_cram), file(row.input_crai)) }
    .set { samples_ch }

samples_ch.into { coverageInChannel; countsInChannel; hsMetricsInChannel; isMetricsInChannel, cnvKitTargetCoverageInChannel; cnvKitAntiTargetCoverageInChannel}

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

process cnvKitTargetCoverage {

    label 'cnvKitTasks'

    input:
    tuple val(sample_id), path(input_cram), path(input_crai) from cnvKitTargetCoverageInChannel

    output:
    tuple val(sample_id), path("*.targetcoverage.cnn") into cnvKitTargetCoverageOutChannel

    script:
    """
    cnvkit.py coverage $input_cram $cnvkit_target_bed -o "$sample_id".targetcoverage.cnn
    """
}

process cnvKitAntiTargetCoverage {

    label 'cnvKitTasks'

    input:
    tuple val(sample_id), path(input_cram), path(input_crai) from cnvKitAntiTargetCoverageInChannel

    output:
    tuple val(sample_id), path("*.antitargetcoverage.cnn") into cnvKitAntiTargetCoverageOutChannel

    script:
    """
    cnvkit.py coverage $input_cram $cnvkit_antitarget_bed  -o "$sample_id".antitargetcoverage.cnn
    """
}


countsOutChannel.into { aggregateCounts_ch; aggregateFpkm_ch }

process aggregateCoverage {
    publishDir "$params.outdir/CombinedCov/", pattern: "*coverage*"

    label 'combineTasks'

    input:
    file input_files from coverageOutChannel.collect()

    output:
    file "${batch}.coverage.bed.gz" into coverage_out_channel

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
    file "${batch}.counts.bed.gz" into counts_out_channel

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
    file "${batch}.fkpm.bed.gz" into fkpm_out_channel

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
    file "${batch}_bioSex-Assignment.txt" into assignedSexChannel

    script:
    """
    assign_bio_sex.py $input_files \
    --suffix ".mosdepth.summary.txt" > "$batch"_bioSex-Assignment.txt
    """
}

process defineProcessGroups {
    publishDir "$params.outdir/ProcessGroups/", pattern: "*.svg"
    publishDir "$params.outdir/ProcessGroups/", pattern: "*stats.txt"

    input:
    file fkpm_file from fkpm_out_channel

    output:
    path "*stats.txt" into defineProcessGroups_out_channel
    path "*100nns.txt" into defineTargetReferences_out_channel
    path "*.svg"

    script:
    """
    define_sub_batches.R $fkpm_file
    """

}

process prepareXhmmInput {
    publishDir "$params.outdir/XHMM_metrics/", pattern: "*XHMM.samples.txt"
    input:
    path coverage_file from coverage_out_channel
    path clusters_file from defineProcessGroups_out_channel

    output:
    file '*xhmm.in.txt' into xhmmInChannel
    file "*XHMM.samples.txt"

    script:
    """
    make_xhmm_input.R $coverage_file $clusters_file
    """
}

process runXhmm {
    publishDir "$params.outdir/XHMM_metrics/", pattern: "*filtered_targets.txt"
    publishDir "$params.outdir/XHMM_metrics/", pattern: "*filtered_samples.txt"
    publishDir "$params.outdir/XHMM_metrics/", pattern: "*PCA.PC.txt"
    publishDir "$params.outdir/XHMM_metrics/", pattern: "*PCA.PC_LOADINGS.txt"
    publishDir "$params.outdir/XHMM_metrics/", pattern: "*PCA.PC_SD.txt"
    publishDir "$params.outdir/XHMM_metrics/", pattern: "*num_removed_PC.txt"
    publishDir "$params.outdir/XHMM_calls/", pattern: "*calls.vcf"
    publishDir "$params.outdir/XHMM_calls/", pattern: "*xcnv"

    input:
    file coverage_file from xhmmInChannel.flatten()

    output:
    path "*filtered_targets.txt"
    path "*filtered_samples.txt"
    path "*.PCA.PC.txt"
    path "*.PCA.PC_LOADINGS.txt"
    path "*.PCA.PC_SD.txt"
    path "*.num_removed_PC.txt"
    path "*calls.vcf"
    path "*xcnv"

    script:
    """
    run_xhmm.sh $coverage_file
    """
    // posteriors?
}


defineTargetReferences_out_channel.into { exomeDepthReferences;  exomeDepthReferencesXchr }


process runExomeDepth{
    publishDir "$params.outdir/ExomeDepth_calls/", pattern: "*calls.bed"
    publishDir "$params.outdir/ExomeDepth_refsets/", pattern: "*.reference.txt"

    input:
    file counts_file from counts_out_channel
    file reference_set from exomeDepthReferences.flatten()

    output:
    path "*calls.bed"
    path "*.reference.txt"

    script:
    """
    ExomeDepth.R $counts_file $reference_set
    """
}


process runExomeDepthXchr{
    publishDir "$params.outdir/ExomeDepth_calls_Xchr/", pattern: "*calls.bed"
    publishDir "$params.outdir/ExomeDepth_refsets_Xchr/", pattern: "*.reference.txt"

    input:
    file bio_sex from assignedSexChannel
    file counts_file from counts_out_channel
    file reference_set from exomeDepthReferencesXchr.flatten()

    output:
    path "*calls.bed"
    path "*.reference.txt"

    script:
    """
    ExomeDepthX.R $counts_file $reference_set $bio_sex
    """
}


/*
Channel
      .fromPath("$projectDir/assets/clamms_special_regions.grch38.bed")
      .set {clammsSpecialRegionsChannel}

process makeCLAMMSWindows {

    input:
    path special_regions from clammsSpecialRegionsChannel

    output:
    path "windows.bed"

    script:
    """
    sed -e 's/^chr//g' $target_bed | awk > target_nochr.bed;


    """

}
*/
