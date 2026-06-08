process CUSTOM_VCF_FILTER {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.8.5--pyhdfd78af_0' :
        'biocontainers/artic:1.8.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}.pass.vcf.gz"), path("${meta.id}.pass.vcf.gz.tbi"), emit: pass_vcf
    tuple val(meta), path("${meta.id}.fail.vcf"), emit: fail_vcf
    path "versions.yml", emit: versions

    script:
    def filterArg = '--nanopolish'
    def argsList = []
    if ( params.variant_caller == "medaka" ) {
        filterArg = "--medaka"
    } else if ( params.variant_caller == "clair3" ) {
        filterArg = "--clair3"
        argsList.add("--min-depth ${params.min_depth}")
        argsList.add("--min-qual-c3 ${params.min_qual_clair3}")
        argsList.add("--min-frameshift-qual ${params.min_frameshift_qual}")
        argsList.add("--min-allele-freq ${params.min_allele_freq}")
        argsList.add("--min-mask-freq ${params.min_mask_freq}")
    }
    def argsConfig = argsList.join(" ")
    """
    cs_vcf_filter.py \\
        $filterArg \\
        $argsConfig \\
        $vcf \\
        ${meta.id}.pass.vcf \\
        ${meta.id}.fail.vcf
    bgzip -f ${meta.id}.pass.vcf
    tabix -p vcf ${meta.id}.pass.vcf.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.pass.vcf.gz
    touch ${meta.id}.pass.vcf.gz.tbi
    touch ${meta.id}.fail.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
