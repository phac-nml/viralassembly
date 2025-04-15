process CHOPPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.7.0--hdcf5f25_0':
        'biocontainers/chopper:0.7.0--hdcf5f25_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.processed.fastq.gz"), emit: fastq
    path "versions.yml", emit: versions

    script:
    // Checking if gzipped or not for stdin to chopper
    //  Note that pipeline should always be just cat
    def cat_cmd = "cat"
    if ( fastq.endsWith('.gz') ) {
        cat_cmd = "zcat"
    }
    """
    $cat_cmd \\
        $fastq \\
    | chopper \\
        --threads $task.cpus \\
        --quality 8 \\
        --minlength 100 \\
    | gzip > ${meta.id}.processed.fastq.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """

    stub:
    """
    # Adding a read to get the option to pass the filtering read check
    #  if we want to
    read="@read1
    TTT
    +
    CCC
    "

    echo -e \$read > ${meta.id}.processed.fastq
    gzip ${meta.id}.processed.fastq

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
