// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_APPLYBQSR {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam), path(bai), path(bqsr_table)
    path(fasta)
    path(fastaidx)
    path(dict)

    output:
    tuple val(meta), path("*.bam"), path("*.bai") , emit: bam
    path "*.version.txt"                          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    gatk ApplyBQSR \\
        -R $fasta \\
        -I $bam \\
        --bqsr-recal-file $bqsr_table \\
        --tmp-dir . \\
        -O ${prefix}.bam \\
        $options.args
    echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > ${software}.version.txt
    """
}
