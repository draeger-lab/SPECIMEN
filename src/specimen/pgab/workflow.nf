#!/usr/bin/env nextflow

params.config_path = 'example_config.yaml'
config_path_ch = Channel.of(params.config_path)

process run_PGAB {
    input:
    path config

    output:
    stdout

    script:
    """
    python3 $projectDir/workflow.py \\
        ${config}
    """
}

workflow {
    run_PGAB(config_path_ch)
}