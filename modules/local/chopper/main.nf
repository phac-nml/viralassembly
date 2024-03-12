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
    sampleName = "$meta.id"
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
    | gzip > ${sampleName}.processed.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
