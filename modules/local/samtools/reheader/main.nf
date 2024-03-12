process SAMTOOLS_REHEADER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    val command

    output:
    tuple val(meta), path("${sampleName}.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    """
    samtools \\
        reheader \\
        ${command} \\
        $bam \\
        > ${sampleName}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
