process ARTIC_ALIGN_TRIM {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.8.5--pyhdfd78af_0' :
        'biocontainers/artic:1.8.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed
    val mode

    output:
    tuple val(meta), path("${meta.id}.*trimmed.rg.sorted.bam"), path("${meta.id}.*trimmed.rg.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    def argsList = []
    if ( params.normalise ) {
        argsList.add("--normalise ${params.normalise}")
    } else {
        argsList.add("--normalise 0")
    }
    outName = "${meta.id}.trimmed.rg.sorted.bam"
    // Start mode = Trim to start of primers instead of ends
    if ( mode == "primers" ) {
        outName = "${meta.id}.primertrimmed.rg.sorted.bam"
    } else {
        argsList.add("--no-trim-primers")
    }
    def argsConfig = argsList.join(" ")
    """
    align_trim \\
        $argsConfig \\
        --report ${meta.id}.alignreport-${mode}.csv \\
        --amp-depth-report ${meta.id}.amplicon_depths.tsv \\
        --primer-match-threshold 15 \\
        $primer_bed \\
        < $bam \\
    | samtools sort -T ${meta.id} - -o $outName

    samtools index $outName

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        align_trim: \$(echo \$(align_trim --version 2>&1) | sed 's/align_trim //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.trimmed.rg.sorted.bam
    touch ${meta.id}.trimmed.rg.sorted.bam.bai

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        align_trim: \$(echo \$(align_trim --version 2>&1) | sed 's/align_trim //')
    END_VERSIONS
    """
}
