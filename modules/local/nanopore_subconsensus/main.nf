/*
    Variant calling at the subconsensus level for nanopore data
        Includes some post processing scripts to reformat/simplify vcfs
*/
process CLAIRSTO_VARIANTS {
    label 'process_medium'
    maxRetries    = 2
    errorStrategy { task.exitStatus == 64 
    ? 'exit' 
    : 'retry' }
    tag "${meta.id}"

    // conda: No functional conda environment currently
    container 'docker://hkubal/clairs-to:v0.4.2'

    input:
    tuple val(meta), path(bam), path(bai)
    path reference
    path fai

    output:
    tuple val(meta),
          path("${meta.id}-clairSTO-out/snv.vcf.gz"), path("${meta.id}-clairSTO-out/snv.vcf.gz.tbi"),
          path("${meta.id}-clairSTO-out/indel.vcf.gz"), path("${meta.id}-clairSTO-out/indel.vcf.gz.tbi"), emit: vcf
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
        --platform ${params.clairS_model} \\
        --output_dir "${meta.id}-clairSTO-out" \\
        -s ${meta.id} \\
        --chunk_size 1000 \\
        --threads 6 \\
        --include_all_ctgs \\
        --disable_verdic \\

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

// combine the indels and snvs into a single VCF
process CAT_VCF {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta),
          path(snv_vcf), path(snv_index),
          path(indel_vcf), path(indel_index)

    output:
    tuple val(meta), path("${meta.id}-cat.vcf.gz"), path("${meta.id}-cat.vcf.gz.tbi"), emit: vcf

    script:
    """
    real_snv=\$(readlink -f $snv_vcf)
    real_indel=\$(readlink -f $indel_vcf)

    # Need to use absolute paths in bcftools command
    bcftools concat \\
        --output-type z \\
        -a -o ${meta.id}-cat.vcf.gz \\
        -O z \$real_snv \$real_indel \\
        
    tabix -p vcf ${meta.id}-cat.vcf.gz
    """
}

// compare with main vcf and remove duplicates
process DEDUP_VCFS {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta),
          path(cat_vcf), path(cat_index),
          path(pass_vcf)

    output:
    // dedup_vcfs/0000.vcf	for records private to	sample-cat.vcf.gz
    tuple val(meta), path("dedup_vcfs/0000.vcf"), emit: vcf

    script:
    """
    real_cat=\$(readlink -f $cat_vcf)
    real_pass=\$(readlink -f $pass_vcf)

    #Need to use absolute paths in bcftools command
    dedup_vcfs.py \\
        --medaka-vcf \$real_pass \\
        --clairs-vcf \$real_cat
    """
}

// fix filters and quality scores
process FIX_VCF {
    label 'process_single'
    tag "$meta.id"
    publishDir "${params.outdir}/vcf", pattern: "*.vcf.gz", mode: "copy"
    publishDir "${params.outdir}/vcf", pattern: "*.vcf.gz.tbi", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(dedup_vcf)

    output:
    tuple val(meta), path("${meta.id}-minorvar.vcf.gz"), path("${meta.id}-minorvar.vcf.gz.tbi"),  emit: vcf

    script:
    """
    fix_clair_vcf.py \\
        -i $dedup_vcf \\
        -o ${meta.id}-minorvar.vcf
    """
}
