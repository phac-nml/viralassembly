process NANOSTAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanostat:1.6.0--pyhdfd78af_0' :
        'biocontainers/nanostat:1.6.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.nanostat.txt"), emit: stats
    path "versions.yml", emit: versions

    script:
    """
    NanoStat \\
        --fastq $fastq \\
        > ${meta.id}.nanostat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NanoStat: \$(NanoStat --version | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.nanostat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        NanoStat: \$(NanoStat --version | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
