process CLAIR3_VARIANTS {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.2.0--py310h779eee5_0' :
        'biocontainers/clair3:1.2.0--py310h779eee5_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path fai
    path model

    output:
    tuple val(meta), path("${meta.id}.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    run_clair3.sh \
        --bam_fn=$bam \\
        --ref_fn=$reference \\
        --threads=${task.cpus} \\
        --platform='ont' \\
        --model_path=$model \\
        --output="${meta.id}-out" \\
        --min_coverage=5 \\
        --haploid_precise \\
        --enable_long_indel \\
        --include_all_ctgs \\
        --chunk_size=10000 \\
        --no_phasing_for_fa

    gunzip ${meta.id}-out/merge_output.vcf.gz
    ln -s ${meta.id}-out/merge_output.vcf ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """
}
