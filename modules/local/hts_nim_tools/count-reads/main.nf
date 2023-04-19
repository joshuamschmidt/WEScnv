process HTSNIMTOOLS_COUNT_READS {
    tag "$meta.id"
    label 'process_medium'

    container 'joshmschmidt/hts_nim_tools:0.2.0'

    input:
    tuple val(meta),  path(bam), path(bai), path(bed)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path('*.cpt.bed.gz')     , emit: counts_bed
    //path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--fasta ${fasta}" : ""
    //def interval = bed ? "--by ${bed}" : ""
    if (bed && args.contains("--by")) {
        exit 1, "'--by' can only be specified once when running mosdepth! Either remove input BED file definition or remove '--by' from 'ext.args' definition"
    }
    if (!bed && args.contains("--thresholds")) {
        exit 1, "'--thresholds' can only be specified in conjunction with '--by'"
    }

    """
    hts_nim_tools count-reads \\
        --threads $task.cpus \\
        $reference \\
        $args \\
        $bed \\
        $bam \\
        | sort -k1,1 -k2,2n \
        | gzip > "$prefix".cpt.bed.gz

    """
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    // END_VERSIONS
    // """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.cpt.bed.gz

    """
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //     mosdepth: \$(mosdepth --version 2>&1 | sed 's/^.*mosdepth //; s/ .*\$//')
    // END_VERSIONS
    // """
}
