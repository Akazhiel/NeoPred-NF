// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_BASERECALIBRATOR_SPARK {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path dict
    path knownSites
    path knownSites_tbi

    output:
    tuple val(meta), path("*.table"), emit: table
    path "*.version.txt" ,            emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def sitesCommand = knownSites.collect{"--known-sites ${it}"}.join(' ')

    """
    gatk BaseRecalibratorSpark  \
        -R $fasta \
        -I $bam \
        $sitesCommand \
        --tmp-dir . \
        $options.args \
        -O ${prefix}.table \
        --spark-master local[${task.cpus}]
    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
}
