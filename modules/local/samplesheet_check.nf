// Import generic module functions
include { saveFiles } from './functions'

params.options = [:]

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'pipeline_info', meta:[:], publish_by_meta:[]) }

    input:
    path samplesheet

    output:
    path '*.csv'

    script: // This script is bundled with the pipeline, in nf-core/neoprednf/bin/
    """
    check_samplesheet.py \\
        $samplesheet \\
        samplesheet.valid.csv
    """
}
