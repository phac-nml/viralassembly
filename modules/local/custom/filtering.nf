// To track filtered out samples
process TRACK_FILTERED_SAMPLES {
    label 'process_single'
    tag "$meta.id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/coreutils:8.31--h14c3975_0'
        : 'biocontainers/coreutils:8.31--h14c3975_0' }"

    input:
    tuple val(meta), path(fastq)
    path metadata
    val fail_message

    output:
    path "*.status.csv", emit: csv

    script:
    """
    SAMPLENAME="${meta.id}"
    IRIDANAME="${meta.irida_id}"

    # Get sample name if metadata is given
    if [ -f "${metadata}" ]; then
        ## Make sure we have barcode and sample column indexes
        barcode_col=\$(awk -v RS='\t' '/barcode/{print NR; exit}' ${metadata})
        sample_col=\$(awk -v RS='\t' '/sample/{print NR; exit}' ${metadata})
        if [ "\$barcode_col" == "" -o "\$sample_col" == "" ]; then
            echo "ERROR: Column 'barcode' or 'sample' does not exist and is required"
            exit 1
        fi

        ## Get sample name based on barcode number - duplicate barcode column should be checked before we get here in pipeline
        barcode_n="\${SAMPLENAME//[!0-9]/}"
        sample_name=\$(awk -F'\t' -v bcol="\$barcode_col" -v scol="\$sample_col" -v barcode_n="\$barcode_n" '{if(\$bcol == barcode_n) print \$scol}' ${metadata})
        ## Set name only if we find a match
        if [ "\$sample_name" != "" ]; then
            SAMPLENAME=\$sample_name
            IRIDANAME=\$sample_name
        fi
    fi

    # Output
    echo "sample,qc_pass,irida_id" > \${SAMPLENAME}.status.csv
    echo "\$SAMPLENAME,${fail_message},\$IRIDANAME" >> \${SAMPLENAME}.status.csv
    """

    stub:
    """
    touch ${meta.id}.status.csv
    """
}
