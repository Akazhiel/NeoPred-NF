// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SOMATICSNIPER {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor)
    path  fasta
    path  fai

    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    bam-somaticsniper \\
        $options.args \\
        -f $fasta \\
        $bam_tumor \\
        $bam_normal \\
        ${prefix}.vcf

    echo '1.0.5.0' > ${software}.version.txt
    """
}
