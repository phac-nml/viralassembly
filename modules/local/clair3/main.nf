process CLAIR3_VARIANTS {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}${pool ? "-${pool}": ""}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:2.0.2--py311hbc58adc_0' :
        'biocontainers/clair3:2.0.2--py311hbc58adc_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(pool), path(pool_bed)
    path reference
    path fai
    path model
    val no_pool_split

    output:
    tuple val(meta), path("${meta.id}${pool ? ".${pool}" : ""}.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''

    // Check for amplicon vs shotgun based variant calling
    def has_pool = pool != null && pool.toString()

    // Define arguments needed if doing amplicon based variant calling
    def bed_args = has_pool ?
        """\
        --bed_fn="${pool_bed}" \\
        --ref_pct_full=1 \\
        --var_pct_full=1 \\
        --enable_variant_calling_at_sequence_head_and_tail
        """
        : ''

    // Specify output name based on amplicon vs shotgun
    def output_vcf = "${meta.id}${has_pool ? ".${pool}" : ""}.vcf"

    // If amplicon, check if pool split bam files is specified
    def split_pool = has_pool && !no_pool_split

    // Based on if were splitting, the bam file has a different name
    def in_bam = split_pool ? "${pool}.sorted.bam" : bam

    // Splitting Command: Filter bam to be the pool only otherwise the full
    //  alignment mode may see the differences in depth at the overlap
    //  region and not call specific SNPs
    def split_command = split_pool ?
        """
        samtools view -b -r "${pool}" "${bam}" -o "${pool}.sorted.bam"
        samtools index "${pool}.sorted.bam"
        """
        : ''

    """
    ${split_command}

    run_clair3.sh \
        $args \\
        --bam_fn=$in_bam \\
        --ref_fn="\$PWD/$reference" \\
        --threads=${task.cpus} \\
        --model_path="$model" \\
        --output="${meta.id}-out" \\
        ${bed_args}

    gunzip ${meta.id}-out/merge_output.vcf.gz
    ln -s ${meta.id}-out/merge_output.vcf "${output_vcf}"

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """

    stub:
    def has_pool = pool != null && pool.toString()
    def output_vcf = "${meta.id}${has_pool ? ".${pool}" : ""}.vcf"
    """
    touch "${output_vcf}"

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """
}
