// Fixes clairS-TO filters and quality scores for viral samples (i.e. removes euk filters and restores QS to QUAL field)
process FIX_VCF {
    label 'process_single'
    label 'error_retry'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(dedup_vcf)

    output:
    tuple val(meta), path("${meta.id}-minor.vcf.gz"),  emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    fix_clair_vcf.py \\
        -i $dedup_vcf \\
        -o ${meta.id}-minor.vcf \\
        -q ${params.min_qual_clairS}

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """
}
