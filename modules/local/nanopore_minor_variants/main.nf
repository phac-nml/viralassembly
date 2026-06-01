/*
    Variant calling at the subconsensus level (i.e. minor variants) from nanopore data
        Includes some post processing scripts to reformat vcfs
*/
process CLAIRSTO_VARIANTS {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}"
    // Only Docker container available - no conda or singularity support
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://hkubal/clairs-to:v0.4.2' :
        'docker.io/hkubal/clairs-to:v0.4.2' }"

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

// combine the indels and snvs into a single VCF
process CAT_VCF {
    label 'process_single'
    label 'error_retry'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(snv_vcf), path(indel_vcf)

    output:
    tuple val(meta), path("${meta.id}-cat.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    real_snv=\$(readlink -f $snv_vcf)
    real_indel=\$(readlink -f $indel_vcf)

    # Need to use absolute paths in bcftools command because links were causing errors.
    bcftools concat \\
        -a -o ${meta.id}-cat.vcf \\
        \$real_snv \$real_indel

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

// compare with main vcf and remove duplicates
process DEDUP_VCFS {
    label 'process_single'
    label 'error_retry'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta),
        path(cat_vcf),
        path(cat_tbi),
        path(pass_vcf)

    output:
    // dedup_vcfs/0000.vcf - records private to sample-cat.vcf.gz (unique to ClairS-TO)
    tuple val(meta), path("dedup_vcfs/0000.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    real_cat=\$(readlink -f $cat_vcf)
    real_pass=\$(readlink -f $pass_vcf)

    #Need to use absolute paths in bcftools command
    dedup_vcfs.py \\
        --consensus-vcf \$real_pass \\
        --clairSTO-vcf \$real_cat

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}

// fix filters and quality scores
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

// combine the major and minor snps/indels into a single VCF
process PASS_VCF {
    label 'process_single'
    label 'error_retry'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(full_vcf)

    output:
    tuple val(meta), path("${meta.id}-full.vcf.gz"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    bcftools view --output-type z -f PASS $full_vcf > ${meta.id}-full.vcf.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
