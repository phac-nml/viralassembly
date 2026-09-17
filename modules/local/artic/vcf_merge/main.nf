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
process ARTIC_VCF_MERGE {
    label 'process_single'
    tag "$meta.id"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.11.1--pyhdfd78af_0' :
        'biocontainers/artic:1.11.1--pyhdfd78af_0' }"

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
