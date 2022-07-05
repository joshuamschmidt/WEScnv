#!/usr/bin/env nextflow

// default params
params.outdir = 'run/'
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
bait_bed=params.bait_bed
cnvkit_target_bed=params.cnvKitTarget
cnvkit_antitarget_bed=params.cnvKitAntiTarget

Channel
    .fromPath(params.inputFile)
    .splitCsv(header:true, sep:'\t')
    .map{ row-> tuple(row.sample_id, file(row.input_cram), file(row.input_crai), file(row.reference_fasta), file(row.reference_idx)) }
    .set { samples_ch }

samples_ch.into { coverageInChannel; countsInChannel; cnvKitTargetCoverageInChannel; cnvKitAntiTargetCoverageInChannel }

process cramCoverage {
    publishDir "$params.outdir/CoverageSummary", pattern: "*.summary.txt"
    publishDir "$params.outdir/Coverage", pattern: "*.regions.bed.gz"

    label 'bamTasks'

    input:
    tuple val(sample_id), path(input_cram), path(input_crai), path(input_reference), path(input_idx) from coverageInChannel

    output:
    path "${sample_id}.regions.bed.gz" into coverageOutChannel
    path "${sample_id}.mosdepth.summary.txt" into coverageSummaryChannel

    script:
    """
    mosdepth --fasta $input_reference \
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
    tuple val(sample_id), path(input_cram), path(input_crai), path(input_reference), path(input_idx) from countsInChannel

    output:
    path "*.cpt.bed.gz" into countsOutChannel

    script:
    """
    hts_nim_tools count-reads \
    --fasta $input_reference \
    --mapq 25 \
    --threads $task.cpus \
    $target_bed $input_cram \
    | sort -k1,1 -k2,2n \
    | gzip > "$sample_id".cpt.bed.gz
    """
}

process cnvKitTargetCoverage {

    input:
    tuple val(sample_id), path(input_cram), path(input_crai), path(input_reference), path(input_idx) from cnvKitTargetCoverageInChannel

    output:
    tuple val(sample_id), path("*.targetcoverage.cnn") into cnvKitTargetCoverageOutChannel

    script:
    """
    mosdepth --fasta $input_reference \
    --by $cnvkit_target_bed \
    --no-per-base \
    --mapq 25 \
    --threads $task.cpus \
    $sample_id \
    $input_cram
    convert_to_cnvkit_coverage.sh "$sample_id".regions.bed.gz $sample_id 'target'
    """
}

process cnvKitAntiTargetCoverage {

    input:
    tuple val(sample_id), path(input_cram), path(input_crai) from cnvKitAntiTargetCoverageInChannel

    output:
    tuple val(sample_id), path("*.antitargetcoverage.cnn") into cnvKitAntiTargetCoverageOutChannel

    script:
    """
    cp $reference_fasta_index .
    mosdepth --fasta $reference_fasta \
    --by $cnvkit_antitarget_bed \
    --no-per-base \
    --mapq 25 \
    --threads $task.cpus \
    $sample_id \
    $input_cram
    convert_to_cnvkit_coverage.sh "$sample_id".regions.bed.gz $sample_id 'anti_target'
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


defineTargetReferences_out_channel.into { exomeDepthReferences;  exomeDepthReferencesXchr; cnvkitReferences }


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


// define all downstream cnvkit channels and processes

cnvKitTargetCoverageOutChannel.into { cnvKitTargetSampleCh; cnvKitTargetRefCh ; cnvKitTargetFixCh }
cnvKitAntiTargetCoverageOutChannel.into { cnvKitAntiTargetSampleCh; cnvKitAntiTargetRefCh ; cnvKitAntiTargetFixCh }

process makeCnvKitSampleRefpairs {

    input:
    file reference_set from cnvkitReferences.flatten()

    output:
    tuple stdout, path(sample_refs) into cnvKitSampleRefCh

    shell:
    '''
    awk 'NR==2 { printf("%s", $1); exit }' !{reference_set}
    awk 'NR!=2 { print $0 }' !{reference_set} > sample_refs
    '''
}

cnvKitTargetSampleCh
    .join(cnvKitAntiTargetSampleCh)
    .set{ cnvKitCombinedSampleCh }

/*
cnvKitCombinedSampleCh
    .join(cnvKitSampleRefCh)
    .set{ makeCnvRefPanelsInCh }

makeCnvRefPanelsInCh
    .first()
    .view()

*/
cnvKitSampleRefCh
    .set{ cnvKitSampleRef_PanelCh }

process makeCnvRefPanels {

    input:
    tuple val(sample_id), path(sample_refs) from cnvKitSampleRef_PanelCh
    file input_target_files from cnvKitTargetRefCh.collect()
    file input_anti_target_files from cnvKitAntiTargetRefCh.collect()

    output:
    tuple val(sample_id), path("*100.cnr") into cnvKitPanelRefCh

    script:
    """
    cp $reference_fasta_index .
    makeCnvKitReference.py --sample_id $sample_id \
    --matched_ref $sample_refs \
    --refFasta $reference_fasta \
    --target_files $input_target_files \
    --anti_target_files $input_anti_target_files
    """
}


/*
process cnvKitFixSample {

    input:
    tuple val(sample_id), path(target_coverage), path(antitarget_coverage), path(reference) from cnvKitPanelRefCh

    output:
    tuple val(sample_id), path("*.cnr") into cnvKitSegmentCh

    script:
    """
    cnvkit.py fix $target_coverage $antitarget_coverage $reference -o "$sample_id".cnr
    """
}
*/
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
