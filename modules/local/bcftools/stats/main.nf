process BCFTOOLS_STATS {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_0' :
        'biocontainers/bcftools:1.19--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf)
    path reference

    output:
    tuple val(meta), path("*stats.txt"), emit: stats
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    """
    # Command #
    bcftools stats \\
        --fasta-ref $reference \\
        $vcf \\
        > ${sampleName}.stats.txt

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
