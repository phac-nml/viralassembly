process PANGOLIN_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pangolin:4.3.4--pyhdfd78af_1' :
        'biocontainers/pangolin:4.3.4--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)
    path(db)

    output:
    tuple val(meta), path('*.csv'), emit: report
    path  "versions.yml"          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def db_command = db ? "--datadir ${db}" : ''
    """
    export XDG_CACHE_HOME=/tmp/.cache

    pangolin \\
        $fasta\\
        $db_command \\
        --outfile ${prefix}.pangolin.csv \\
        --threads $task.cpus \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangolin: \$(pangolin --version | sed "s/pangolin //g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export XDG_CACHE_HOME=/tmp/.cache

    touch ${prefix}.pangolin.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pangolin: \$(pangolin --version | sed "s/pangolin //g")
    END_VERSIONS
    """
}
