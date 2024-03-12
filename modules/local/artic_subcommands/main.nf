/*
    Artic subcommands and some custom adaptations of them to work with clair3 without
        explicitly having to use the minion.py workflow
*/

// Helper function for combine VCFs in the format needed for artic merge
def transformVCFList (inputList) {
    def transformedOutput = inputList.collect { entry ->
        "${entry[1]}:${entry[0]}"
    }.join(" ")
    return transformedOutput
}

// Subcommands start here
process ARTIC_ALIGN_TRIM {
    label 'process_single'
    tag "$meta.id"
    publishDir "${params.outdir}/bam", pattern: "${sampleName}.*trimmed.rg.sorted.bam", mode: "copy"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${sampleName}.alignreport-*", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path primer_bed
    val mode

    output:
    tuple val(meta), path("${sampleName}.*trimmed.rg.sorted.bam"), path("${sampleName}.*trimmed.rg.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    def argsList = []
    if ( params.normalise ) {
        argsList.add("--normalise ${params.normalise}")
    }
    outName = "${sampleName}.primertrimmed.rg.sorted.bam"
    // Start mode = Trim to start of primers instead of ends
    if ( mode == "start" ) {
        outName = "${sampleName}.trimmed.rg.sorted.bam"
        argsList.add("--start")
    }
    def argsConfig = argsList.join(" ")
    """
    align_trim \\
        $argsConfig \\
        --remove-incorrect-pairs \\
        --report ${sampleName}.alignreport-${mode}.txt \\
        $primer_bed \\
        < $bam 2> ${sampleName}.alignreport-${mode}.er | samtools sort -T ${sampleName} - -o $outName

    samtools index $outName

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
process ARTIC_VCF_MERGE {
    label 'process_single'
    tag "$meta.id"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${sampleName}.merged.vcf", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    // The vcf_tuples input is [[ path(vcf), val(pool) ], [...]]
    //   The path(vcf) is turned into a string of the full path using the val() input type
    //   The process still works, just is a bit iffy I'd say
    input:
    tuple val(meta), val(vcf_tuples)
    path primer_bed

    output:
    tuple val(meta), path("${sampleName}.merged.vcf"), emit: vcf
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    def vcfs = transformVCFList(vcf_tuples)
    """
    artic_vcf_merge \\
        ${sampleName} \\
        $primer_bed \\
        2> ${sampleName}.primersitereport.txt \\
        $vcfs

    # Versions from nf-core #
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
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${sampleName}*.vcf.gz"), path("${sampleName}*.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    """
    bgzip -f $vcf
    tabix -f -p vcf ${vcf}.gz

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
process CUSTOM_VCF_FILTER {
    label 'process_single'
    tag "$meta.id"
    publishDir "${params.outdir}/vcf", pattern: "${sampleName}.pass.vcf.gz*", mode: "copy"
    publishDir "${params.outdir}/vcf", pattern: "${sampleName}.fail.vcf", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("${sampleName}.pass.vcf.gz"), path("${sampleName}.pass.vcf.gz.tbi"), emit: pass_vcf
    tuple val(meta), path("${sampleName}.fail.vcf"), emit: fail_vcf
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    def filterArg = '--nanopolish'
    if ( params.variant_caller == 'medaka' ) {
        filterArg = '--medaka'
    } else if ( params.variant_caller == 'clair3' ) {
        filterArg = '--clair3'
    }
    """
    cs_vcf_filter.py \\
        $filterArg \\
        $vcf \\
        ${sampleName}.pass.vcf \\
        ${sampleName}.fail.vcf
    bgzip -f ${sampleName}.pass.vcf
    tabix -p vcf ${sampleName}.pass.vcf.gz

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
process ARTIC_MAKE_DEPTH_MASK{
    label 'process_single'
    tag "$meta.id"
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${sampleName}.coverage_mask.txt", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path reference

    output:
    tuple val(meta), path("${sampleName}.coverage_mask.txt"), emit: coverage_mask
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    """
    artic_make_depth_mask \\
        --store-rg-depths \\
        $reference \\
        $bam \\
        ${sampleName}.coverage_mask.txt

    # Versions from nf-core #
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
    // publishDir "${params.outdir}/articMinionNextflow", pattern: "${sampleName}.coverage_mask.txt", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path(reference)

    output:
    tuple val(meta), path("${sampleName}.coverage_mask.txt"), emit: coverage_mask

    script:
    sampleName = "$meta.id"
    """
    cs_make_depth_mask.py \\
        --depth 20 \\
        $reference \\
        $bam \\
        ${sampleName}.coverage_mask.txt
    """
}
process ARTIC_MASK {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(coverage_mask), path(fail_vcf)
    path reference

    output:
    tuple val(meta), path("${sampleName}.preconsensus.fasta"), emit: preconsensus
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    """
    artic_mask \\
        $reference \\
        $coverage_mask \\
        $fail_vcf \\
        ${sampleName}.preconsensus.fasta

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
