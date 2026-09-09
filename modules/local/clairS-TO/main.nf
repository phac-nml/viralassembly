/*
    Variant calling at the subconsensus level (i.e. minor variants) from nanopore data
*/
process CLAIRSTO_VARIANTS {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}"
    // Only Docker container available - no conda or singularity support
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clairs-to:v0.4.4' :
        'docker.io/hkubal/clairs-to:v0.4.4' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path fai
    val model

    output:
    tuple val(meta),
        path("${meta.id}-clairSTO-out/snv.vcf.gz"),
        path("${meta.id}-clairSTO-out/indel.vcf.gz"),
        emit: vcf
    path "versions.yml", emit: versions

    script:

    def argsList = []

    argsList.add("--snv_min_af ${params.min_snv_af_clairS}")
    argsList.add("--indel_min_af ${params.min_indel_af_clairS}")
    argsList.add("--min_coverage ${params.min_cov_clairS}")
    argsList.add("--qual ${params.min_qual_clairS}")

    def argsConfig = argsList.join(" ")

    """
    run_clairs_to \\
        ${argsConfig} \\
        --tumor_bam_fn $bam \\
        --ref_fn $reference \\
        --threads ${task.cpus} \\
        --platform "$model" \\
        --output_dir "${meta.id}-clairSTO-out" \\
        --chunk_size 1000 \\
        --include_all_ctgs \\
        --disable_verdict

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairSTO: \$(echo \$(run_clairs_to -v) | sed 's/run_clairs_to //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        clairSTO: \$(echo \$(run_clairs_to -v) | sed 's/run_clairs_to //')
    END_VERSIONS
    """
}
