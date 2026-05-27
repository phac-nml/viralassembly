process CLAIR3_VARIANTS {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}-${pool}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/clair3:1.2.0--py310h779eee5_0' :
        'biocontainers/clair3:1.2.0--py310h779eee5_0' }"

    input:
    tuple val(meta), path(bam), path(bai), val(pool), path(pool_bed)
    path reference
    path fai
    path model
    val no_pool_split

    output:
    tuple val(meta), path("${meta.id}.${pool}.vcf"), val(pool), emit: vcf
    path "versions.yml", emit: versions

    script:
    // Based on if were splitting the bam file has a different name
    def in_bam = "${pool}.sorted.bam"
    if ( no_pool_split ) {
        in_bam = bam
    }
    """
    # Filter bam to be the pool only when running with pools
    #  Otherwise the full alignment mode may see the differences in depth at the overlap
    #  region and not call specific SNPs
    if [[ "$no_pool_split" == "false" ]]; then
        samtools view -b -r $pool $bam -o ${pool}.sorted.bam
        samtools index ${pool}.sorted.bam
    fi

    run_clair3.sh \
        --bam_fn=$in_bam \\
        --bed_fn=$pool_bed \\
        --ref_fn=$reference \\
        --threads=${task.cpus} \\
        --platform='ont' \\
        --model_path="$model" \\
        --output="${meta.id}-out" \\
        --min_coverage=5 \\
        --haploid_precise \\
        --enable_long_indel \\
        --include_all_ctgs \\
        --ref_pct_full=1 \\
        --var_pct_full=1 \\
        --chunk_size=5000 \\
        --no_phasing_for_fa \\
        --enable_variant_calling_at_sequence_head_and_tail

    gunzip ${meta.id}-out/merge_output.vcf.gz
    ln -s ${meta.id}-out/merge_output.vcf ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.${pool}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clair3: \$(echo \$(run_clair3.sh -v) | sed 's/Clair3 //')
    END_VERSIONS
    """
}
