// QC Modules
//  Using artic env for the moment as it has a bunch of tools and
//  is used in earlier steps already
process MAKE_SAMPLE_QC_CSV {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.11.1--pyhdfd78af_0' :
        'biocontainers/artic:1.11.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(consensus), path(bam), path(bai), path(depth_bed), path(vcf), path(min_vcf), path(nextclade_csv)
    path primer_bed
    path metadata
    path pcr_primers
    val skip_nextclade

    output:
    tuple val(meta), path ("${meta.id}.qc.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    // Need to structure args based on what we have
    def metadataArg = metadata ? "--metadata $metadata" : ""
    def seqArg = primer_bed ? "--seq_bed $primer_bed" : ""
    def pcrArg = pcr_primers ? "--pcr_bed $pcr_primers" : ""
    def minVcfArg = min_vcf ? "--min_vcf $min_vcf" : ""
    def nextcladeColumnsArgs = skip_nextclade ? "" : "--add_nextclade_columns"
    def nextcladeCsvArg = nextclade_csv ? "--nextclade_csv $nextclade_csv" : ""
    """
    qc.py \\
        --bam $bam \\
        --vcf $vcf \\
        --consensus $consensus \\
        --depth $depth_bed \\
        $metadataArg \\
        $seqArg \\
        $pcrArg \\
        $minVcfArg \\
        $nextcladeColumnsArgs \\
        $nextcladeCsvArg \\
        --sample $meta.id \\
        --irida_id $meta.irida_id

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.qc.csv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process FINAL_QC_CSV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.11.1--pyhdfd78af_0' :
        'biocontainers/artic:1.11.1--pyhdfd78af_0' }"

    input:
    path combined_csv
    path filter_tracking_csv
    path metadata
    path reference
    val neg_control_threshold
    val neg_ctrl_substrings

    output:
    path "overall.qc.csv", emit: csv
    path "versions.yml", emit: versions

    script:
    // Need to structure args based on what we have
    def version = workflow.manifest.version
    def filterArg = filter_tracking_csv ? "--filter_tracking $filter_tracking_csv" : ""
    def metadataArg = metadata ? "--metadata $metadata" : ""
    """
    final_checks.py \\
        $filterArg \\
        $metadataArg \\
        --threshold $neg_control_threshold \\
        --neg_ctrl_substrings '$neg_ctrl_substrings' \\
        --csv $combined_csv \\
        --reference $reference \\
        --version $version

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch overall.qc.csv

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
