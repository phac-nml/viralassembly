process LONGSHOT {
    label 'process_medium'
    tag "$meta.id"
    publishDir "${params.outdir}/vcf", pattern: "${meta.id}.longshot.merged.vcf", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longshot:1.0.0--hd4f2111_2':
        'biocontainers/longshot:1.0.0--hd4f2111_2' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(bam), path(bai)
    path reference
    path reference_fai

    output:
    tuple val(meta), path("${meta.id}.longshot.merged.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def VERSION = '1.0.0' // Longshot version does not seem to be being printed out
    """
    longshot \\
        -P 0 \\
        -F \\
        -A \\
        --no_haps \\
        --bam $bam \\
        --ref $reference \\
        --out ${meta.id}.longshot.merged.vcf \\
        --potential_variants $vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longshot: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '1.0.0' // Longshot version does not seem to be being printed out
    """
    touch ${meta.id}.longshot.merged.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longshot: $VERSION
    END_VERSIONS
    """
}
