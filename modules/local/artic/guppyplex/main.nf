process ARTIC_GUPPYPLEX {
    label 'process_medium'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.2--pyhdfd78af_0' :
        'biocontainers/artic:1.6.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.fastq"), emit: fastq
    path("versions.yml"), emit: versions

    script:
    // Fastq input can either be a directory or a set of fastq files
    //  Outputs are the same then after allowing a streamlined pipeline
    if ( fastq.isDirectory() ) {
        """
        artic guppyplex \\
            --min-length ${params.min_length} \\
            --max-length ${params.max_length} \\
            --output ${meta.id}.fastq \\
            --directory $fastq

        # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        END_VERSIONS
        """
    } else {
        """
        mkdir -p input_fastq
        mv $fastq input_fastq/
        artic guppyplex \\
            --min-length ${params.min_length} \\
            --max-length ${params.max_length} \\
            --output ${meta.id}.fastq \\
            --directory input_fastq

        # Versions #
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
        END_VERSIONS
        """
    }

    stub:
    """
    touch ${meta.id}.fastq

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
