// QC Modules
//  Using artic env for the moment as it has a bunch of tools and
//  is used in earlier steps
process MAKE_SAMPLE_QC_CSV {
    label 'process_single'
    tag "$meta.id"
    publishDir "${params.outdir}/sample_csvs", pattern: "*.csv", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    tuple val(meta), path(consensus), path(bam), path(bai), path(depth_bed), path(vcf)
    path primer_bed
    path metadata
    path pcr_primers

    output:
    tuple val(meta), path ("${sampleName}.qc.csv"), emit: csv
    path "versions.yml", emit: versions

    script:
    sampleName = "$meta.id"
    // Need to structure args based on what we have
    def version = workflow.manifest.version
    def metadataArg = metadata ? "--metadata $metadata" : ""
    def seqArg = primer_bed ? "--seq_bed $primer_bed" : ""
    def pcrArg = pcr_primers ? "--pcr_bed $pcr_primers" : ""
    def analysis = "nanopolish"
    if ( params.variant_caller == 'medaka' ) {
        analysis = "medaka"
    } else if ( params.variant_caller == 'clair3' ) {
        analysis = "clair3"
    }
    """
    qc.py \\
        --analysis $analysis \\
        --consensus $consensus \\
        --bam $bam \\
        --vcf $vcf \\
        --depth $depth_bed \\
        $metadataArg \\
        $seqArg \\
        $pcrArg \\
        --sample $sampleName

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
process FINAL_QC_CSV {
    label 'process_single'
    publishDir "${params.outdir}", pattern: "overall.qc.csv", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.4--pyh7cba7a3_1' :
        'biocontainers/artic:1.2.4--pyh7cba7a3_1' }"

    input:
    path combined_csv
    path filter_tracking_csv
    path metadata
    path reference
    val neg_control_threshold
    val neg_ctrl_substrings

    output:
    path "overall.qc.csv", emit: csv
    path "versions.yml", emit: versions

    script:
    // Need to structure args based on what we have
    def version = workflow.manifest.version
    def filterArg = filter_tracking_csv ? "--filter_tracking $filter_tracking_csv" : ""
    def metadataArg = metadata ? "--metadata $metadata" : ""
    """
    final_checks.py \\
        $filterArg \\
        $metadataArg \\
        --threshold $neg_control_threshold \\
        --neg_ctrl_substrings '$neg_ctrl_substrings' \\
        --csv $combined_csv \\
        --reference $reference \\
        --version $version

    # Versions from nf-core #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
