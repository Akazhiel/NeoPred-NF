// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process GATK4_BASERECALIBRATOR {
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
    path knownSites
    path knownSites_tbi
    path target_bed

    output:
    tuple val(meta), path("*.table"), emit: table
    path "*.version.txt" ,            emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def sitesCommand = knownSites.collect{"--known-sites ${it}"}.join(' ')
    def intervals = target_bed ? "--intervals ${target_bed}" : ""

    if (meta.type == "dna") {
    """
    gatk BaseRecalibrator  \
        -R $fasta \
        -I $bam \
        $sitesCommand \
        --tmp-dir . \
        $options.args \
        $intervals \
        -O ${prefix}.table
    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
    } else {
    """
    gatk BaseRecalibrator  \
        -R $fasta \
        -I $bam \
        $sitesCommand \
        --tmp-dir . \
        --use-original-qualities \
        $options.args \
        -O ${prefix}.table
    gatk --version | grep Picard | sed "s/Picard Version: //g" > ${software}.version.txt
    """
    }
    
}
