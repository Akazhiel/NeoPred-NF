// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process YARA_MAPPER {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:[]) }

    input:
    tuple val(meta), path(bam), path(bai)
    path  index

    output:
    tuple val(meta), path("*.mapped.bam") , emit: bam
    path "*.version.txt"                  , emit: version

    script:
    def software = getSoftwareName(task.process)
    def split_cpus = Math.floor(task.cpus/2)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    samtools view -@ ${split_cpus} -h -f 0x40 $bam > output_1.bam
    samtools view -@ ${split_cpus} -h -f 0x80 $bam > output_2.bam
    samtools bam2fq output_1.bam > output_1.fastq
    samtools bam2fq output_2.bam > output_2.fastq
    yara_mapper \\
        $options.args \\
        -t ${task.cpus} \\
        -f bam \\
        ${index}/yara \\
        output_1.fastq \\
        output_2.fastq > output.bam
    samtools view -@ ${split_cpus} -hF 4 -f 0x40 -b output.bam > ${prefix}_1.mapped.bam
    samtools view -@ ${split_cpus} -hF 4 -f 0x80 -b output.bam > ${prefix}_2.mapped.bam
    echo \$(yara_mapper --version) | sed "s/^.*yara_mapper version: //; s/ .*\$//" > ${software}.version.txt
    """
}
