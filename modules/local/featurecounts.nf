// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process FEATURECOUNTS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bams)
    path(annotation)

    output:
    tuple val(meta), path("*featureCounts.txt")        , emit: counts
    tuple val(meta), path("*featureCounts.txt.summary"), emit: summary
    path "*.version.txt"                               , emit: version

    script:
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def paired_end = '-p'

    """
    featureCounts \\
        $options.args \\
        $paired_end \\
        -T $task.cpus \\
        -a $annotation \\
        -o ${prefix}.featureCounts.txt \\
        ${bam}

    echo \$(featureCounts -v 2>&1) | sed 's/featureCounts v//g'> ${software}.version.txt
    """
}
