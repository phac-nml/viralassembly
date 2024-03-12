process SAMTOOLS_DEPTH {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${sampleName}.depth.bed"), emit: bed
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    """
    echo -e "chrom\tpos\tdepth" \\
        > ${sampleName}.depth.bed
    samtools depth \\
        -a \\
        $bam \\
        >> ${sampleName}.depth.bed

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
