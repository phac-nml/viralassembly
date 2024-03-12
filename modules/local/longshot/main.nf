process LONGSHOT {
    label 'process_medium'
    tag "$meta.id"
    publishDir "${params.outdir}/vcf", pattern: "${sampleName}.longshot.merged.vcf", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longshot:0.4.5--hd175d40_2':
        'biocontainers/longshot:0.4.5--hd175d40_2' }"

    input:
    tuple val(meta), path(vcf), path(tbi), path(bam), path(bai)
    path reference
    path reference_fai

    output:
    tuple val(meta), path("${sampleName}.longshot.merged.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    def VERSION = '0.4.5' // Longshot version does not seem to be being printed out
    """
    longshot \\
        -P 0 \\
        -F \\
        -A \\
        --no_haps \\
        --bam $bam \\
        --ref $reference \\
        --out ${sampleName}.longshot.merged.vcf \\
        --potential_variants $vcf

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longshot: $VERSION
    END_VERSIONS
    """
}
