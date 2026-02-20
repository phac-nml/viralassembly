process MEDAKA_CONSENSUS {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.3--py39h05d5c5e_0' :
        'biocontainers/medaka:1.11.3--py39h05d5c5e_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.hdf"), emit: hdf
    path "versions.yml", emit: versions

    script:
    """
    medaka consensus \\
        --model ${params.medaka_model} \\
        --threads ${task.cpus} \\
        --chunk_len 800 \\
        --chunk_ovlp 400 \\
        $bam \\
        ${meta.id}.hdf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.hdf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
process MEDAKA_VARIANT {
    label 'process_medium'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.11.3--py39h05d5c5e_0' :
        'biocontainers/medaka:1.11.3--py39h05d5c5e_0' }"

    input:
    tuple val(meta), path(hdf)
    path reference

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    medaka variant \\
        $reference \\
        $hdf \\
        ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
