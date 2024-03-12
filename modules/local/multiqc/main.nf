/*
    Pipeline specific multiqc processes for sample level and run level reports
*/
process MULTIQC_SAMPLE {
    label 'process_single'
    publishDir "${params.outdir}/sample_mqc", pattern: "*.html", mode: "copy"
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.20--pyhdfd78af_0' :
        'biocontainers/multiqc:1.20--pyhdfd78af_0' }"

    input:
    path multiqc_config
    tuple val(meta), path(sample_csv), path(variation_csv), path(consensus_variant_tsv), path(qualimap_bamqc_data), path(amp_depth_tsv)

    output:
    path "*.html", emit: html

    script:
    sampleName = "$meta.id"
    """
    multiqc \\
        -f \\
        -k yaml \\
        --config $multiqc_config \\
        --filename ${sampleName}.report.html \\
        --title "NML ${sampleName} Sample Report" \\
        .
    """
}
process MULTIQC_OVERALL {
    label 'process_single'
    publishDir "${params.outdir}", pattern: "*.html", mode: "copy"

    conda "${moduleDir}/environment.yml"
    container "https://depot.galaxyproject.org/singularity/multiqc:1.20--pyhdfd78af_2"

    input:
    path multiqc_config
    path amp_depth_tsvs
    path merged_amplicon_comp_csv
    path bcftools_stats
    path samtools_flagstats
    path qualimap_bamqc_data
    path snpeff_csvs
    path qc_csv
    path versions_yml

    output:
    path "*.html", emit: html

    script:
    """
    multiqc -f -k yaml --config $multiqc_config .
    """
}
