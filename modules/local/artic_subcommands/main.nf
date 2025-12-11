/*
    Artic subcommands and some custom adaptations of them to work with clair3 without
        explicitly having to use the minion.py workflow
*/

// Helper function for combine VCFs in the format needed for artic merge
//  Files are always NAME.POOL.vcf based on previous process
//  So if the name has a . the size() - 2 hopefully gets the right number
def transformVCFList (inputList) {
    def transformedOutput = inputList.collect { vcf ->
        def name_split = vcf.name.split(/\./)
        def pool = name_split[name_split.size()-2]
        "${pool}:${vcf}"
    }.join(" ")
    return transformedOutput
}

//
// Subcommands start here
process ARTIC_ALIGN_TRIM {
    label 'process_single'
    tag "$meta.id"
    publishDir "${params.outdir}/bam", pattern: "${meta.id}.*trimmed.rg.sorted.bam*", mode: "copy"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${meta.id}.alignreport-*", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed
    val mode

    output:
    tuple val(meta), path("${meta.id}.*trimmed.rg.sorted.bam"), path("${meta.id}.*trimmed.rg.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    def argsList = []
    if ( params.normalise ) {
        argsList.add("--normalise ${params.normalise}")
    }
    outName = "${meta.id}.trimmed.rg.sorted.bam"
    // Start mode = Trim to start of primers instead of ends
    if ( mode == "primers" ) {
        outName = "${meta.id}.primertrimmed.rg.sorted.bam"
        argsList.add("--trim-primers")
    }
    def argsConfig = argsList.join(" ")
    """
    align_trim \\
        $argsConfig \\
        --remove-incorrect-pairs \\
        --report ${meta.id}.alignreport-${mode}.csv \\
        $primer_bed \\
        < $bam 2> ${meta.id}.alignreport-${mode}.er \\
    | samtools sort -T ${meta.id} - -o $outName

    samtools index $outName

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.trimmed.rg.sorted.bam
    touch ${meta.id}.trimmed.rg.sorted.bam.bai

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process ARTIC_VCF_MERGE {
    label 'process_single'
    tag "$meta.id"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${meta.id}.merged.vcf", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    // The vcf_tuples input is [[ path(vcf), val(pool) ], [...]]
    //   The path(vcf) is turned into a string of the full path using the val() input type
    //   The process still works, just is a bit iffy I'd say
    input:
    tuple val(meta), path(vcfs)
    path primer_bed

    output:
    tuple val(meta), path("${meta.id}.merged.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    def vcfs_in_str = transformVCFList(vcfs)
    """
    artic_vcf_merge \\
        ${meta.id} \\
        $primer_bed \\
        2> ${meta.id}.primersitereport.txt \\
        $vcfs_in_str

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.merged.vcf

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process ZIP_AND_INDEX_VCF {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${meta.id}*.vcf.gz"), path("${meta.id}*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    """
    bgzip -f $vcf
    tabix -f -p vcf ${vcf}.gz

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.vcf.gz
    touch ${meta.id}.vcf.gz.tbi

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process CUSTOM_VCF_FILTER {
    label 'process_single'
    tag "$meta.id"
    publishDir "${params.outdir}/vcf", pattern: "${meta.id}.pass.vcf.gz*", mode: "copy"
    publishDir "${params.outdir}/vcf", pattern: "${meta.id}.fail.vcf", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

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
    }
    """
    cs_vcf_filter.py \\
        $filterArg \\
        $argsList \\
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

process ARTIC_MAKE_DEPTH_MASK{
    label 'process_single'
    tag "$meta.id"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${meta.id}.coverage_mask.txt", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("${meta.id}.coverage_mask.txt"), emit: coverage_mask
    path "versions.yml", emit: versions

    script:
    """
    artic_make_depth_mask \\
        --store-rg-depths \\
        $reference \\
        $bam \\
        ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

// Slow but the bedtools adaptation I was working on I couldn't quite get to be genomic index
//  Will have to look at that more as it was a lot quicker
process CUSTOM_MAKE_DEPTH_MASK {
    label 'process_medium'
    label 'error_retry'
    tag "${meta.id}"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${meta.id}.coverage_mask.txt", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)

    output:
    tuple val(meta), path("${meta.id}.coverage_mask.txt"), emit: coverage_mask

    script:
    """
    cs_make_depth_mask.py \\
        --depth 20 \\
        $reference \\
        $bam \\
        ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.coverage_mask.txt

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}

process ARTIC_MASK {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.7.4--pyhdfd78af_0' :
        'biocontainers/artic:1.7.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(coverage_mask), path(fail_vcf)
    path reference

    output:
    tuple val(meta), path("${meta.id}.preconsensus.fasta"), emit: preconsensus
    path "versions.yml", emit: versions

    script:
    """
    artic_mask \\
        $reference \\
        $coverage_mask \\
        $fail_vcf \\
        ${meta.id}.preconsensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.preconsensus.fasta

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
