// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STAR_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(reads)
    path  index
    path  gtf

    output:
    tuple val(meta), path('*d.out.bam')       , emit: bam
    path "*versions.txt"                      , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def seq_center = "--outSAMattrRGline ID:$prefix PU:$prefix SM:$prefix LB:$prefix PL:ILLUMINA"

    """
    STAR \\
        --genomeDir $index \\
        --readFilesIn $reads  \\
        --runThreadN $task.cpus \\
        --outFileNamePrefix $prefix. \\
        $seq_center \\
        $options.args

    echo \$(STAR --version 2>&1) > ${software}.versions.txt
    """
}
