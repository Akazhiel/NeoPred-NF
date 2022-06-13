// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_MARKDUPLICATES {
    tag "$meta.id"
    label 'process_medium'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam)
    val use_metrics

    output:
    tuple val(meta), path("*.bam"), path("*.bai")      , emit: bam
    tuple val(meta), path("*.metrics"), optional : true, emit: metrics
    path "*.version.txt"                               , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def metrics  = use_metrics ? "M=${prefix}.metrics" :''
    // def bams     = bam.collect(){ x -> "INPUT=".concat(x.toString()) }.join(" ")

    def markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    """
    gatk --java-options ${markdup_java_options} \\
        MarkDuplicates \\
        $metrics \\
        INPUT=$bam \\
        TMP_DIR=. \\
        ASSUME_SORT_ORDER=coordinate \\
        O=${prefix}.bam \\
        $options.args

    mv ${prefix}.bai ${prefix}.bam.bai

    echo \$(gatk MarkDuplicates --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}
