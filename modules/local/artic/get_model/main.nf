process ARTIC_GET_MODELS {
    label 'process_single'
    tag "Download"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.11.1--pyhdfd78af_0' :
        'biocontainers/artic:1.11.1--pyhdfd78af_0' }"

    input:
    val model

    output:
    path("clair3_models/${model}"), emit: model
    path("versions.yml"), emit: versions

    script:
    """
    CLAIR_BIN_DIR=\$(dirname \$(which run_clair3.sh))

    mkdir -p clair3_models
    if ls \$CLAIR_BIN_DIR/models/ | grep -q "$model" ; then
        cp -r \$CLAIR_BIN_DIR/models/$model ./clair3_models/
    else
        artic_get_models --model-dir ./clair3_models --models $model
    fi

    # Confirm pt file exists
    if ! ls ./clair3_models/$model | grep -q "\\.pt" ; then
        echo "ERROR: Missing .pt file in downloaded $model folder. Check model name, artic version, and download source to try to resolve"
        exit 1
    fi

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p clair3_models/$model

    # Versions #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(echo \$(artic --version 2>&1) | sed 's/artic //')
    END_VERSIONS
    """
}
