process PRIMALBEDTOOLS_AMPLICON {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/primalbedtools:1.0.0--pyhdfd78af_0' :
        'biocontainers/primalbedtools:1.0.0--pyhdfd78af_0' }"

    input:
    path bed

    output:
    path "amplicon.bed", emit: bed
    path "versions.yml", emit: versions

    script:
    """
    primalbedtools \\
        amplicon \\
        $bed \\
        --primertrim \\
    > amplicon.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Primalbedtools: \$(primalbedtools --version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    touch amplicon.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Primalbedtools: \$(primalbedtools --version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
