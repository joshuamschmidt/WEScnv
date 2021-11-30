#!/usr/bin/env nextflow

// default params
params.usePicard = false
params.useGATK3 = false
params.aggMetrics = false
params.outputDir = "run/"
params.saveMode = "copy"
params.type = "WGS"
params.cohort = "general"
params.target_intervals = false
params.bait_intervals = false


inputFile = file(params.inputFile)
lines = inputFile.readLines()
cohort = params.cohort
type = params.type
manifest = params.manifest

bait_intervals = params.bait_intervals
usePicard = params.usePicard
useGATK3 = params.useGATK3
aggMetrics = params.aggMetrics
saveMode = params.saveMode
target_intervals = params.target_intervals

usrgitc = "/usr/gitc/"


process mosdepthCounts {

    publishDir "${params.outputDir}/CoverageSummary/", mode: "${saveMode}" ,pattern: "*mosdepth.summary.txt"

    input:
    each line from lines

    output:
    set val(ius), val(processing_chain), file("${output_bam_basename}.bam"), file("${output_bam_basename}.bam.bai") into FastqAndBwaAndSorMadupOutChannel
    set val(ius), file ("${metrics_filename}") into MarkDupMetrics
    file "${metrics_filename}" into MultiQCDupsInput

    script:
    list = line.split('\t')
    ius = list[0]
    input_R1 = list[1]
    input_R2 = list[2]
    INPUT_RGID = list[3]
    INPUT_RGPL = list[4]
    INPUT_RGPU = list[5]
    INPUT_RGLB = list[6]
    INPUT_RGSM = list[7]
    INPUT_RGCN = list[8]
    processing_chain = "aligned.sorted.markduped"
    output_bam_basename = "${ius}.${processing_chain}"
    metrics_filename = "${ius}.duplicate_metrics"
    """
    cp $input_R1 \${TMPDIR}/
    cp $input_R2 \${TMPDIR}/

    bwa mem \
        -t 16 \
        -R @RG\\\\tID:"${INPUT_RGID}"\\\\tPL:"${INPUT_RGPL}"\\\\tPU:"${INPUT_RGPU}"\\\\tLB:"${INPUT_RGLB}"\\\\tSM:"${INPUT_RGSM}" \
        -K 100000000 \
        -v 3 -Y \$ref_fasta \
        \${TMPDIR}/\$(basename ${input_R1}) \${TMPDIR}/\$(basename ${input_R2}) | \
    bamsormadup \
        threads=16 \
        inputformat=sam \
        outputformat=bam \
        reference=\${ref_fasta} \
        M=${metrics_filename} \
        indexfilename="${output_bam_basename}.bam.bai" \
        optminpixeldif=2500 \
        O=${output_bam_basename}.bam > ${output_bam_basename}.bam
    """
}

FastqAndBwaAndSorMadupOutChannel.set { SorMaDupChannel }


process CreateSequenceGroupingTSV {

    output:
    file "sequence_grouping.txt" into SeqGroupChannel
    file "sequence_grouping_with_unmapped.txt" into SeqGroupUnMChannel

    script:
    """
#!/usr/bin/env python3
with open("/directflow/ClinicalGenomicsPipeline/dev/2019-11-05-BINF-469-GATK4/gatk-workflows/broad-references/hg38/v0/Homo_sapiens_assembly38.dict", "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\\t")
            # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
# the last element after a :, so we add this as a sacrificial element.
hg38_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += "\\t" + sequence_tuple[0] + hg38_protection_tag
    else:
        tsv_string += "\\n" + sequence_tuple[0] + hg38_protection_tag
        temp_size = sequence_tuple[1]
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
with open("sequence_grouping.txt","w") as tsv_file:
    tsv_file.write(tsv_string)
    tsv_file.close()

tsv_string += '\\n' + "unmapped"
with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
    tsv_file_with_unmapped.write(tsv_string)
    tsv_file_with_unmapped.close()
    """
}

SorMaDupChannel.into { CheckContInChannel ; BaseRecalibratorInChannel ; BQSRInChannel }


if(target_intervals){

    process CheckContaminationWES {

        label 'Contamination'

        publishDir "${params.outputDir}/verifyBAMID/", mode: "${saveMode}"

        input:
        set val(ius), val(processing_chain), file(input_bam), file(input_bai) from CheckContInChannel
        set file(exome_contamination_sites_ud), file(exome_contamination_sites_bed), file(exome_contamination_sites_mu) from SubsettedContamResources

        output:
        set val(ius), file("${output_prefix}.selfSM") into getValInChannelWES

        script:
        output_prefix = "${ius}.preBqsr"
        """
        set -e

        # creates a ${output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
        # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
        /usr/gitc/VerifyBamID \
            --Verbose \
            --NumPC 4 \
            --Output ${output_prefix} \
            --BamFile ${input_bam} \
            --Reference \${ref_fasta} \
            --UDPath ${exome_contamination_sites_ud} \
            --MeanPath ${exome_contamination_sites_mu} \
            --BedPath ${exome_contamination_sites_bed} \
            1>/dev/null
        """
    }

    getValInChannelWES.set { getValInChannel }

} else {
    process CheckContamination {

        label 'Contamination'

        publishDir "${params.outputDir}/verifyBAMID/", mode: "${saveMode}"

        input:
        set val(ius), val(processing_chain), file(input_bam), file(input_bai) from CheckContInChannel

        output:
        set val(ius), file("${output_prefix}.selfSM") into getValInChannelWGS

        script:
        output_prefix = "${ius}.preBqsr"
        """
        set -e

        # creates a ${output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
        # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
        /usr/gitc/VerifyBamID \
            --Verbose \
            --NumPC 4 \
            --Output ${output_prefix} \
            --BamFile ${input_bam} \
            --Reference \${ref_fasta} \
            --UDPath \${contamination_sites_ud} \
            --MeanPath \${contamination_sites_mu} \
            --BedPath \${contamination_sites_bed} \
        1>/dev/null
        """
    }

    getValInChannelWGS.set { getValInChannel }

}

getValInChannel.into { getValInChannel ; renameBamIn }

process CheckContaminationPython {

    input:
    set val(ius), file(selfSM) from getValInChannel

    output:
    set val(ius), stdout into CheckContaminationOut

    script:
    contamination_underestimation_factor = 0.75
    """
#!/usr/bin/env python3
import csv
import sys
with open('${selfSM}') as selfSM:
    reader = csv.DictReader(selfSM, delimiter='\t')
    i = 0
    for row in reader:
        if i != 0:
            if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
                # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
                # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
                # vcf and bam.
                sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
                sys.exit(1)
            print(float(row["FREEMIX"])/0.75)
        i = i + 1
    # there should be exactly two rows, and if this isn't the case the format of the output is unexpectedly different
    # and the results are not reliable.
    if i != 2:
        sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
        sys.exit(2)
    """
}

CheckContaminationOut.into { CheckContaminationOut ; ContaminationFileIn }

SeqGroupChannel.set { RecalSG }

RecalRange = Channel.from( 0..17 )

process BaseRecalibrator {

    label 'QC'

    input:
    set val(ius), val(processing_chain), file(input_bam), file(input_bai) from BaseRecalibratorInChannel
    each i from RecalRange
    val sequence_group_interval_all from RecalSG.readLines()

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

BaseRecalOut.groupTuple(size:18).set { GatherBQSRInChannel }


process GatherBQSRReports {

    publishDir "${params.outputDir}/BQSRReport/", mode: "${saveMode}"

    input:
    set val(ius), file(input_bqsr_reports) from GatherBQSRInChannel

    output:
    set val(ius), file("${output_report_filename}") into GatherBQSROutChannel
    file("${output_report_filename}") into MultiQCBQSRInput

    script:
    output_report_filename = "${ius}.recal_data.csv"
    """
    gatk --java-options "-Xms3000m" \
        GatherBQSRReports \
        -I \$(echo $input_bqsr_reports | sed 's/ /\\n/g' | sort | tr '\\n' ' ' | sed 's/ / -I /g' | sed 's/ -I \$//' ) \
        -O ${output_report_filename}
    """

}

SeqGroupUnMChannel.set { BQSRSeqGUnMChannel }

BQSRInChannel.join(GatherBQSROutChannel, remainder:true).set { BQSRInChannel }

BQSRRange = Channel.from( 0..18 )

process ApplyBQSR {

    input:
    set val(ius), val(processing_chain), file(input_bam), file(input_bai), file(recalibration_report) from BQSRInChannel
    val sequence_group_interval_all from BQSRSeqGUnMChannel.readLines()
    each i from BQSRRange

    output:
    set val(ius), file("${output_bam_basename}.bam"), file("${output_bam_basename}.bai") into BQSROutchannel
    val processing_chain into FinalProcessingChain

    script:
    processing_chain = "${processing_chain}.recalibrated"
    sequence_group_interval = sequence_group_interval_all.get(i).toString()
    if ( i < 10 ) {
        output_bam_basename = "${ius}.${processing_chain}.0${i}"
    } else {
        output_bam_basename = "${ius}.${processing_chain}.${i}"
    }

    """
    gatk --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps \
        -XX:+PrintGCDetails -Xloggc:gc_log.log \
        -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Dsamjdk.compression_level=\${compression_level} -Xms3000m" \
        ApplyBQSR \
        --create-output-bam-md5 \
        --add-output-sam-program-record \
        -R \${ref_fasta} \
        -I ${input_bam} \
        --use-original-qualities \
        -O ${output_bam_basename}.bam \
        -bqsr ${recalibration_report} \
        --static-quantized-quals 10 \
        --static-quantized-quals 20 \
        --static-quantized-quals 30 \
        -L \$(echo ${sequence_group_interval} | sed 's/ / -L /g')
    """
}

BQSROutchannel.groupTuple(size:19).set { GatherBamsInChannel }


process GatherBamFiles {

    publishDir "${params.outputDir}/inputBam/", mode: "symlink", pattern: "*bam*"

    input:
    set val(ius), file(input_bams), file(input_bais) from GatherBamsInChannel

    output:
    set val(ius), file("${ius}.bam"), file("${ius}.bai") into GatherBamsOutChannel

    script:
    """
    java -Dsamjdk.compression_level=\${compression_level} -Xms2000m -jar ${usrgitc}/picard.jar \
        GatherBamFiles \
        INPUT=\$(echo "${input_bams}" | sed 's/ /\\n/g' | sort | tr '\\n' ' ' | sed 's/ / INPUT= /g' | sed 's/ INPUT= \$//' ) \
        OUTPUT=${ius}.bam \
        CREATE_INDEX=true \
        CREATE_MD5_FILE=true
    """
}


GatherBamsOutChannel.into { CramInputChannel ; HCBamsInputChannel ; BamQCMetricsChannel  }

process ConvertToCram {

    publishDir "${params.outputDir}/inputCram/", mode: "symlink", pattern: "*cram*"

    input:
    set val(ius), file(input_bam), file(input_bai) from CramInputChannel

    output:
    set val(ius), file("${output_basename}.cram"), file("${output_basename}.cram.crai") into CramOutputChannel
    set val(ius), file("${output_basename}.cram"), file("${output_basename}.cram.crai"), file("${output_basename}.cram.md5") into GCPCramInput
    file "${output_basename}.cram.md5"

    script:
    output_basename = "${ius}"
    """
    set -e
    set -o pipefail

    samtools view -C -T \${ref_fasta} ${input_bam} | \
    tee ${output_basename}.cram | \
    md5sum | awk '{print \$1}' > ${output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    seq_cache_populate.pl -root ./ref/cache \${ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ${output_basename}.cram
    """
}

// CheckContaminationOut.println()
// GatherBamsOutChannel.println()

HCBamsInputChannel.join(CheckContaminationOut, remainder:true).set { HCInChannel }

ScatterRange = Channel.from( 1..50 )

// HCInChannel.println()
if (useGATK3) {

    process HaplotypeCallerGATK3 {

        input:
        set val(ius), file(input_bam), file(input_bai), val(contamination) from HCInChannel
        each i from ScatterRange

        output:
        set val(ius), file("${ius}*.vcf.gz"), file("${ius}*.vcf.gz.tbi") into HCOutChannel

        script:
        """
        interval_list=\$(head -n${i} \${scatter_order} | tail -n1)
        number=\$(echo \$(basename \$(dirname \$interval_list)) | cut -d_ -f2)
        gvcf_basename="${ius}.\${number}.g"
        contami=\$(echo "${contamination}" | tr -d '\\040\\011\\012\\015')
        ${usrgitc}/gatk4/gatk-launch --javaOptions "-Xms2g" \
            PrintReads \
            -I ${input_bam} \
            --interval_padding 500 \
            -L ${interval_list} \
            -O local.sharded.bam \
         && \
        java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m \
            -jar ${usrgitc}/GATK35.jar \
            -T HaplotypeCaller \
            -R \${ref_fasta} \
            -o \${gvcf_basename}.vcf.gz \
            -I local.sharded.bam \
            -L ${interval_list} \
            -ERC GVCF \
            --max_alternate_alleles 3 \
            -variant_index_parameter 128000 \
            -variant_index_type LINEAR \
            -contamination \${contami} \
            --read_filter OverclippedRead
        """
    }

} else {

    process HaplotypeCaller {
        input:
        set val(ius), file(input_bam), file(input_bai), val(contamination) from HCInChannel
        each i from ScatterRange

        output:
        set val(ius), file("${ius}*.vcf.gz"), file("${ius}*.vcf.gz.tbi") into HCOutChannel

        script:
        output_file_name = ""
        """
        set -e
        interval_list=\$TMPDIR/\$(head -n${i} \${scatter_order} | tail -n1)
        number=\$(echo \$(basename \$(dirname \$interval_list)) | cut -d_ -f2)
        gvcf_basename="${ius}.\${number}.g.vcf.gz"
        contami=\$(echo "${contamination}" | tr -d '\\040\\011\\012\\015')
        gatk --java-options "-Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" \
            HaplotypeCaller \
            -R \${ref_fasta} \
            -I ${input_bam} \
            -L \${interval_list} \
            -O \${gvcf_basename} \
            -contamination \${contami} \
            -G StandardAnnotation -G StandardHCAnnotation -G AS_StandardAnnotation \
            -new-qual \
            -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
            -ERC GVCF
        """
    }

}

HCOutChannel.groupTuple(size:50).set { MergeVCFSInChannel }

process MergeVCFs {

    label 'QC'

    publishDir "${params.outputDir}/inputGvcf/", mode: "${saveMode}", pattern: "*vcf.gz*"

    input:
    set val(ius), file(input_gvcfs), file(input_tbis) from MergeVCFSInChannel

    output:
    set val(ius), file("${output_vcf_name}"), file("${output_vcf_name}.tbi") into MergeVCFSOutChannel

    script:
    output_vcf_name = "${ius}.g.vcf.gz"
    """
    java -Xms2000m -jar ${usrgitc}/picard.jar \
        MergeVcfs \
        INPUT=\$(echo "${input_gvcfs}" | sed 's/ /\\n/g' | sort | tr '\\n' ' ' | sed 's/ / INPUT= /g' | sed 's/ INPUT= \$//' ) \
        OUTPUT=${output_vcf_name}
    """
}


if(target_intervals) {

    BamQCMetricsChannel.into { CollectHsMetricsInputChannel ; CollectRGQCInputChannel ; CollectAggregationMetricsInput ; CollectWgsMetricsInput ; CollectRawWgsMetricsInput ; CalculateRGChecksumInput }

    process CollectHsMetrics {
        label 'QC'

        publishDir "${params.outputDir}/CollectHsMetrics/", mode: "${saveMode}"

        input:
        set val(ius), file(input_bam), file(input_bai) from CollectHsMetricsInputChannel
        file target_intervals

        output:
        file "${metrics_filename}*"

        script:
        metrics_filename = "${ius}.CollectHsMetrics"
        """
        java -Xms5000m -jar ${usrgitc}/picard.jar \
            CollectHsMetrics \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=\${ref_fasta} \
            VALIDATION_STRINGENCY=SILENT \
            TARGET_INTERVALS=${target_intervals} \
            BAIT_INTERVALS=${bait_intervals} \
            METRIC_ACCUMULATION_LEVEL=null \
            METRIC_ACCUMULATION_LEVEL=SAMPLE \
            METRIC_ACCUMULATION_LEVEL=LIBRARY \
            OUTPUT=${metrics_filename}
        """
    }

} else {
    BamQCMetricsChannel.into { CollectRGQCInputChannel ; CollectAggregationMetricsInput ; CollectWgsMetricsInput ; CollectRawWgsMetricsInput ; CalculateRGChecksumInput }
}

process CollectReadgroupBamQualityMetrics {

    label 'QC'

    publishDir "${params.outputDir}/CollectRGBamQC/", mode: "${saveMode}"

    input:
    set val(ius), file(input_bam), file(input_bai) from CollectRGQCInputChannel

    output:
    set file("${output_bam_prefix}.alignment_summary_metrics"), file("${output_bam_prefix}.gc_bias.detail_metrics"), file("${output_bam_prefix}.gc_bias.pdf"), file("${output_bam_prefix}.gc_bias.summary_metrics") into MultiQCCollectRGQCInput

    script:
    output_bam_prefix = "${ius}.readgroup"
    """
    touch ${output_bam_prefix}.gc_bias.detail_metrics \
        ${output_bam_prefix}.gc_bias.pdf \
        ${output_bam_prefix}.gc_bias.summary_metrics

    java -Xms5000m -jar ${usrgitc}/picard.jar \
        CollectMultipleMetrics \
        INPUT=${input_bam} \
        REFERENCE_SEQUENCE=\${ref_fasta} \
        OUTPUT=${output_bam_prefix} \
        ASSUME_SORTED=true \
        PROGRAM="null" \
        PROGRAM="CollectAlignmentSummaryMetrics" \
        PROGRAM="CollectGcBiasMetrics" \
        METRIC_ACCUMULATION_LEVEL="null" \
        METRIC_ACCUMULATION_LEVEL="READ_GROUP"
    """
}

if (aggMetrics) {
    process CollectAggregationMetrics {
        label 'QC'

        publishDir "${params.outputDir}/CollectAggregationMetrics/", mode: "${saveMode}"

        input:
        set val(ius), file(input_bam), file(input_bai) from CollectAggregationMetricsInput

        output:
        set val(ius), file("${output_bam_prefix}.alignment_summary_metrics") into CollectAggSummary
        file("${output_bam_prefix}.alignment_summary_metrics") into MultiQCAggInput
        set file("${output_bam_prefix}.bait_bias_detail_metrics"), file("${output_bam_prefix}.bait_bias_summary_metrics"), file("${output_bam_prefix}.pre_adapter_detail_metrics"), file("${output_bam_prefix}.pre_adapter_summary_metrics") into MultiQCAlignmentInput
        set file("${output_bam_prefix}.gc_bias.detail_metrics"), file("${output_bam_prefix}.gc_bias.pdf"), file("${output_bam_prefix}.gc_bias.summary_metrics") into MultiQCGCbiasInupt
        set file("${output_bam_prefix}.insert_size_histogram.pdf"), file("${output_bam_prefix}.insert_size_metrics") into MultiQCInsertSizeInput
        set file("${output_bam_prefix}.quality_distribution.pdf"), file("${output_bam_prefix}.quality_distribution_metrics") into MultiQCQualInput

        script:
        output_bam_prefix = "${ius}"
        """
        touch ${output_bam_prefix}.gc_bias.detail_metrics \
            ${output_bam_prefix}.gc_bias.pdf \
            ${output_bam_prefix}.gc_bias.summary_metrics \
            ${output_bam_prefix}.insert_size_metrics \
            ${output_bam_prefix}.insert_size_histogram.pdf


        java -Xms5000m -jar ${usrgitc}/picard.jar \
            CollectMultipleMetrics \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=\${ref_fasta} \
            OUTPUT=${output_bam_prefix} \
            ASSUME_SORTED=true \
            PROGRAM="null" \
            PROGRAM="CollectAlignmentSummaryMetrics" \
            PROGRAM="CollectInsertSizeMetrics" \
            PROGRAM="CollectSequencingArtifactMetrics" \
            PROGRAM="CollectGcBiasMetrics" \
            PROGRAM="QualityScoreDistribution" \
            METRIC_ACCUMULATION_LEVEL="null" \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY"
        """
    }

} else {
    CollectAggregationMetricsInput.into { CollectAlignmentSummaryMetricsInput ; CollectInsertSizeMetricsInput ; CollectSequencingArtifactMetricsInput ; CollectGcBiasMetricsInput ;  QualityScoreDistributionInput }

    process CollectAlignmentSummaryMetrics {

        label 'QC'

        publishDir "${params.outputDir}/CollectAlignmentSummaryMetrics/", mode: "${saveMode}"

        input:
        set val(ius), file(input_bam), file(input_bai) from CollectAlignmentSummaryMetricsInput

        output:
        set val(ius), file("${output_bam_prefix}") into CollectAggSummary
        file("${output_bam_prefix}") into MultiQCAggInput

        script:
        output_bam_prefix = "${ius}.alignment_summary_metrics"
        """
        java -Xms1000m -jar ${usrgitc}/picard.jar \
            CollectAlignmentSummaryMetrics \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=\${ref_fasta} \
            OUTPUT=${output_bam_prefix} \
            ASSUME_SORTED=true \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY"
        """
    }

    process CollectInsertSizeMetrics {

        label 'QC'

        publishDir "${params.outputDir}/CollectInsertSizeMetrics/", mode: "${saveMode}"

        input:
        set val(ius), file(input_bam), file(input_bai) from CollectInsertSizeMetricsInput

        output:
        set file("${output_bam_prefix}_histogram.pdf"), file("${output_bam_prefix}_metrics") into MultiQCInsertSizeInput

        script:
        output_bam_prefix = "${ius}.insert_size"
        """
        java -Xms1000m -jar ${usrgitc}/picard.jar \
            CollectInsertSizeMetrics \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=\${ref_fasta} \
            OUTPUT=${output_bam_prefix}_metrics \
            ASSUME_SORTED=true \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY" \
            HISTOGRAM_FILE=${output_bam_prefix}_histogram.pdf

        """
    }

    process CollectSequencingArtifactMetrics {

        label 'QC'

        publishDir "${params.outputDir}/CollectSequencingArtifactMetrics/", mode: "${saveMode}"

        input:
        set val(ius), file(input_bam), file(input_bai) from CollectSequencingArtifactMetricsInput

        output:
        set file("${output_bam_prefix}.bait_bias_detail_metrics"), file("${output_bam_prefix}.bait_bias_summary_metrics"), file("${output_bam_prefix}.pre_adapter_detail_metrics"), file("${output_bam_prefix}.pre_adapter_summary_metrics") into MultiQCAlignmentInput

        script:
        output_bam_prefix = "${ius}"
        """
        java -Xms1000m -jar ${usrgitc}/picard.jar \
            CollectSequencingArtifactMetrics \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=\${ref_fasta} \
            OUTPUT=${output_bam_prefix} \
            ASSUME_SORTED=true
        """
    }

    process CollectGcBiasMetrics {

        label 'QC'

        publishDir "${params.outputDir}/CollectGcBiasMetrics/", mode: "${saveMode}"

        input:
        set val(ius), file(input_bam), file(input_bai) from CollectGcBiasMetricsInput

        output:
        set file("${output_bam_prefix}.detail_metrics"), file("${output_bam_prefix}.pdf"), file("${output_bam_prefix}.summary_metrics") into MultiQCGCbiasInupt

        script:
        output_bam_prefix = "${ius}.gc_bias"
        """
        java -Xms1000m -jar ${usrgitc}/picard.jar \
            CollectGcBiasMetrics \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=\${ref_fasta} \
            OUTPUT=${output_bam_prefix}.detail_metrics \
            ASSUME_SORTED=true \
            METRIC_ACCUMULATION_LEVEL="SAMPLE" \
            METRIC_ACCUMULATION_LEVEL="LIBRARY" \
            CHART_OUTPUT=${output_bam_prefix}.pdf \
            SUMMARY_OUTPUT=${output_bam_prefix}.summary_metrics

        """
    }

    process QualityScoreDistribution {

        label 'QC'

        publishDir "${params.outputDir}/QualityScoreDistribution/", mode: "${saveMode}"

        input:
        set val(ius), file(input_bam), file(input_bai) from QualityScoreDistributionInput

        output:
        set file("${output_bam_prefix}.pdf"), file("${output_bam_prefix}_metrics") into MultiQCQualInput

        script:
        output_bam_prefix = "${ius}.quality_distribution"
        """
        java -Xms1000m -jar ${usrgitc}/picard.jar \
            QualityScoreDistribution \
            INPUT=${input_bam} \
            REFERENCE_SEQUENCE=\${ref_fasta} \
            OUTPUT=${output_bam_prefix}_metrics \
            ASSUME_SORTED=true \
            CHART_OUTPUT=${output_bam_prefix}.pdf
        """
    }
}

process CollectWgsMetrics {

    label 'QC'

    publishDir "${params.outputDir}/CollectWgsMetrics/", mode: "${saveMode}"

    input:
    set val(ius), file(input_bam), file(input_bai) from CollectWgsMetricsInput

    output:
    file("${metrics_filename}") into MultiQCWGSInput

    script:
    metrics_filename = "${ius}.wgs_metrics"
    """
    java -Xms2000m -jar ${usrgitc}/picard.jar \
        CollectWgsMetrics \
        INPUT=${input_bam} \
        VALIDATION_STRINGENCY=SILENT \
        REFERENCE_SEQUENCE=\${ref_fasta} \
        INCLUDE_BQ_HISTOGRAM=true \
        INTERVALS=\${wgs_coverage_interval_list} \
        OUTPUT=${metrics_filename} \
        USE_FAST_ALGORITHM=true \
        READ_LENGTH=250
    """
}

process CollectRawWgsMetrics {

    label 'QC'

    publishDir "${params.outputDir}/CollectRawWgsMetrics/", mode: "${saveMode}"

    input:
    set val(ius), file(input_bam), file(input_bai) from CollectRawWgsMetricsInput

    output:
    file("${metrics_filename}") into MultiQCRawWGSInput

    script:
    metrics_filename = "${ius}.raw_wgs_metrics"
    """
    java -Xms2000m -jar ${usrgitc}/picard.jar \
        CollectRawWgsMetrics \
        INPUT=${input_bam} \
        VALIDATION_STRINGENCY=SILENT \
        REFERENCE_SEQUENCE=\${ref_fasta} \
        INCLUDE_BQ_HISTOGRAM=true \
        INTERVALS=${wgs_coverage_interval_list} \
        OUTPUT=${metrics_filename} \
        USE_FAST_ALGORITHM=true \
        READ_LENGTH=250
    """
}


MarkDupMetrics.join(CollectAggSummary, remainder:true).set { CheckPreValidationInput }

process CheckPreValidation {

    input:
    set val(ius), file(duplication_metrics), file(chimerism_metrics) from CheckPreValidationInput

    output:
    set val(ius), file('duplication.csv'), file('chimerism.csv') into CheckPreValOut

    script:
    """
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ${duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ${chimerism_metrics} | grep -v OF_PAIR > chimerism.csv
    """
}

process CheckPreValidationPython {

    input:
    set val(ius), file(duplication), file(chimerism) from CheckPreValOut

    output:
    set val(ius), stdout into IsOutlierDataOut

    script:
    max_duplication_in_reasonable_sample = 0.30
    max_chimerism_in_reasonable_sample = 0.15
    """
#!/usr/bin/env python3

import csv
with open('duplication.csv') as dupfile:
    reader = csv.DictReader(dupfile, delimiter='\t')
    for row in reader:
        with open("duplication_value.txt","w") as file:
            file.write(row['PERCENT_DUPLICATION'])
            file.close()

with open('chimerism.csv') as chimfile:
    reader = csv.DictReader(chimfile, delimiter='\t')
    for row in reader:
        with open("chimerism_value.txt","w") as file:
            file.write(row['PCT_CHIMERAS'])
            file.close()

# custom code
dups = open("duplication_value.txt", "r")
dups_content = dups.read()
duplication_rate = float(dups_content)

chim = open("chimerism_value.txt", "r")
chim_content = chim.read()
chimerism_rate = float(chim_content)

if (duplication_rate > ${max_duplication_in_reasonable_sample} or chimerism_rate > ${max_chimerism_in_reasonable_sample}):
    print("true")
else:
    print("false")
    """
}

CramOutputChannel.join(IsOutlierDataOut, remainder:true).set { ValidateSamInput }

if (usePicard) {
    process ValidateSamFile {

        label 'QC'

        input:
        set val(ius), file(input_bam), file(input_bai), val(is_outlier_data) from ValidateSamInput

        output:
        file "${report_filename}"

        script:
        report_filename = "${ius}.cram.validation_report"

        """
        java -Xms6000m -jar ${usrgitc}/picard.jar \
            ValidateSamFile \
            INPUT=${input_bam} \
            OUTPUT=${report_filename} \
            REFERENCE_SEQUENCE=\${ref_fasta} \
            MAX_OUTPUT=1000000000 \
            IGNORE="MISSING_TAG_NM" \
            MODE=VERBOSE \
            SKIP_MATE_VALIDATION=\$(echo "${is_outlier_data}" | tr -d '\\040\\011\\012\\015') \
            IS_BISULFITE_SEQUENCED=false
        """
    }
} else {
    process bamvalidate {

        input:
        set val(ius), file(input_bam), file(input_bai), val(is_outlier_data) from ValidateSamInput

        output:
        file "${report_filename}"

        script:
        report_filename = "${ius}.cram.validation_report"
        """
        touch "${report_filename}"

        bamvalidate \
            inputformat=cram \
            I=${input_bam} \
            reference=\${ref_fasta} \
            basequalhist=1 \
            O=WhatIsThis.QuestionMark
        """
    }

}


MergeVCFSOutChannel.into { ValidateGVCFInput ; CollectGvcfCallingMetricsInput ; GCPGvcfInput }

process ValidateGVCF {

    input:
    set val(ius), file(input_vcf), file(input_tbi) from ValidateGVCFInput

    script:
    """
    gatk --java-options "-Xms3000m" \
        ValidateVariants \
        -V ${input_vcf} \
        -R \${ref_fasta} \
        -L \${wgs_calling_interval_list} \
        -gvcf \
        --validation-type-to-exclude ALLELES \
        --dbsnp \${dbSNP_vcf}
    """
}



process CollectGvcfCallingMetrics {

    label 'QC'

    publishDir "${params.outputDir}/CollectGvcfCallingMetrics/", mode: "${saveMode}"

    input:
    set val(ius), file(input_vcf), file(input_tbi) from CollectGvcfCallingMetricsInput

    output:
    set file("${metrics_basename}.variant_calling_summary_metrics"), file("${metrics_basename}.variant_calling_detail_metrics") into MultiQCCollectgvcfInput

    script:
    metrics_basename = "${ius}.g.vcf"
    """
    java -Xms2000m -jar ${usrgitc}/picard.jar \
        CollectVariantCallingMetrics \
        INPUT=${input_vcf} \
        OUTPUT=${metrics_basename} \
        DBSNP=\${dbSNP_vcf} \
        SEQUENCE_DICTIONARY=\${ref_dict} \
        TARGET_INTERVALS=\${wgs_evaluation_interval_list} \
        GVCF_INPUT=true
    """
}



process verifyBamIDMultiQCEdit {

    input:
    set val(ius), val(contam) from ContaminationFileIn

    output:
    file "${ius}.selfSM" into FalconQCContamiInput

    script:
    """
    contamination=\$(echo "${contam}" | tr -d '\\040\\011\\012\\015')
    awk -vp=\$contamination -vq=0.75 'BEGIN{printf "%.10f" ,p * q}' > ${ius}.selfSM
    """

}


// MuitQC and FalconQC
process MultiQC {
    publishDir "${params.outputDir}/multiqc/", mode: "${saveMode}" // , pattern: "*html"

    input:
    file variant_calling_summary_metrics from MultiQCCollectgvcfInput.collect()
    file raw from MultiQCRawWGSInput.collect()
    file wgs from MultiQCWGSInput.collect()
    file bait_bias_detail_metrics from MultiQCAlignmentInput.collect()
    file gc_bias_detail_metrics from MultiQCGCbiasInupt.collect()
    file insert_size_histogram_pdf from MultiQCInsertSizeInput.collect()
    file quality_distribution_pdf from MultiQCQualInput.collect()
    file alignment_summary_metrics from MultiQCCollectRGQCInput.collect()
    file bqsr from MultiQCBQSRInput.collect()
    file dups from MultiQCDupsInput.collect()

    output:
    set file("multiqc_report.html"), file("multiqc_data") into FalconQCInput

    script:
    """
    echo -e -n "biobambam2/bamsormadup:\\n   contents: '# bamsormadup'\\n   \\n   num_lines: 2\\nfn_clean_exts:\\n   - '.duplicate_metrics'\\n   - '.recal_data'\\n   - '.bam'\\n   - '.g.vcf.gz'" > multiqc.config
    multiqc ./ -c multiqc.config
    """
}

process FlaconQC {

    publishDir "${params.outputDir}/FalconQC/", mode: "${saveMode}"

    input:
    file inputFile
    set file(report), file(data) from FalconQCInput
    file contamination from FalconQCContamiInput.collect()

    output:
    file "${manifest}.csv" into GCPUploadQCInput

    script:
    """
    echo "Sample Name,Cohort Name,Batch Name,Flowcell.Lane,Library ID,Platform,Centre of Sequencing,Reference Genome,Type,Description" > sample_metadata.csv
    for sample in \$(cut -d\$'\t' -f1 $inputFile ); do
        platform=\$(grep \$sample $inputFile | cut -d\$'\\t' -f5 | sed 's/\\t/,/g')
        FCIDlane=\$(grep \$sample $inputFile | cut -d\$'\\t' -f6 | sed 's/\\t/,/g')
        lib=\$(grep \$sample $inputFile | cut -d\$'\\t' -f7 | sed 's/\\t/,/g')
        centre=\$(grep \$sample $inputFile | cut -d\$'\\t' -f9 | sed 's/\\t/,/g')
        contamination=\$(cat \${sample}.selfSM)
        echo "\$sample,$cohort,$manifest,\$FCIDlane,\$lib,\$platform,\$centre,hg38,$type,\\"FREEMIX:\$contamination\\""
    done >> sample_metadata.csv

    falcon_multiqc save --directory ./ --sample_metadata sample_metadata.csv --batch_description $manifest

    falcon_multiqc query --select sample --select tool-metric --cohort $cohort \
        --tool-metric picard_AlignmentSummaryMetrics PCT_CHIMERAS 0 0 \
        --tool-metric dup_bamsormadup PERCENT_DUPLICATION 0 0 \
        --tool-metric picard_insertSize MEDIAN_INSERT_SIZE 0 0 \
        --tool-metric picard_wgsmetrics MEDIAN_COVERAGE 0 0 \
        --output ./ --csv --filename ${manifest}

    sed -i 's/sample.description/raw_data.FREEMIX/' ${manifest}.csv
    sed -i 's/FREEMIX://' ${manifest}.csv
    """
}
