// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWA_MEM {
    tag "$meta.id"
    label 'process_high'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(reads)
    path  index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "*.version.txt"         , emit: version

    script:
    def split_cpus = Math.floor(task.cpus/2)
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}.${part}" : "${meta.id}."
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
    
    println "Mapping ${meta.patient}"

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`
    bwa mem \\
        -t ${split_cpus} \\
        $options.args \\
        $read_group \\
        \$INDEX \\
        $reads \\
        | samtools $options.args2 --threads ${split_cpus} -o ${prefix}bam
    echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//' > ${software}.version.txt
    """
}
