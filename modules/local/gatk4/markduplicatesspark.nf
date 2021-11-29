// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options    = initOptions(params.options)

process GATK4_MARKDUPLICATES_SPARK {
    tag "$meta.id"
    label 'process_high'

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam)
    path(reference)
    path(dict) //need to be present in the path
    path(fai)  //need to be present in the path

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: output
    path("*.version.txt")                              , emit: version

    script:
    def software = getSoftwareName(task.process)
    // def bams = bam.collect(){ x -> "-I ".concat(x.toString()) }.join(" ")
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def markdup_java_options = (task.memory.toGiga() > 8) ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicatesSpark \
        -I $bam \
        -O ${prefix}.bam \
        --reference ${reference} \
        --tmp-dir . \
        --spark-master local[${task.cpus}] \\
        $options.args
    echo \$(gatk MarkDuplicatesSpark --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}
