// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process GET_SOFTWARE_VERSIONS {
    publishDir "${params.outdir}/${meta.patient}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:meta, publish_by_meta: false) }

    container "quay.io/biocontainers/python:3.8.3"

    cache false

    input:
    tuple val(meta), path(reads)
    path versions

    output:
    path "software_versions.tsv"     , emit: tsv
    path 'software_versions_mqc.yaml', emit: yaml

    script: // This script is bundled with the pipeline, in nf-core/neoprednf/bin/
    """
    echo $workflow.manifest.version > pipeline.version.txt
    echo $workflow.nextflow.version > nextflow.version.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}
