// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process SPLITNCIGAR {
    tag "$meta.id"
    label 'process_high'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict

    output:
    tuple val(meta), path("*.bam"), path("*.bai"), emit: bam
    path "*.version.txt"                         , emit: versions

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    gatk SplitNCigarReads  \
        --create-output-bam-index \
        -R $fasta \
        -I $bam \
        -O ${prefix}.bam
        
    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
}