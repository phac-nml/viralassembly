process BCFTOOLS_NORM {
    label 'process_single'
    tag "$meta.id"
    publishDir "${params.outdir}/vcf", pattern: "${sampleName}.pass.norm.vcf.gz", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.19--h8b25389_0' :
        'biocontainers/bcftools:1.19--h8b25389_0' }"

    input:
    tuple val(meta), path(preconsensus), path(pass_vcf), path(pass_vcf_tbi)

    output:
    tuple val(meta), path("${sampleName}.pass.norm.vcf.gz"), path("${sampleName}.pass.norm.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    """
    # Fixes variants that are in both the pass and fail vcf that were masked #
    bcftools norm \\
        --check-ref s \\
        -f $preconsensus \\
        $pass_vcf \\
        > ${sampleName}.pass.norm.vcf
    bgzip ${sampleName}.pass.norm.vcf
    tabix ${sampleName}.pass.norm.vcf.gz

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
