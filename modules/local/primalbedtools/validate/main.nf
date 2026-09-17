process PRIMALBEDTOOLS_VALIDATE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/primalbedtools:1.0.0--pyhdfd78af_0' :
        'biocontainers/primalbedtools:1.0.0--pyhdfd78af_0' }"

    input:
    path bed
    path reference

    output:
    path "versions.yml", emit: versions

    script:
    """
    # Validate bed formatting
    primalbedtools \\
        validate_bedfile \\
        $bed

    # Validate bed to ref
    primalbedtools \\
        validate \\
        $bed \\
        $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Primalbedtools: \$(primalbedtools --version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Primalbedtools: \$(primalbedtools --version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
