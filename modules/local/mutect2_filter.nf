// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MUTECT2_FILTER {
    tag "Mutect_filtering"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(mutect_vcf)
    tuple val(meta), path(mutect_vcf_stats)
    path fasta
    path fasta_fai
    path dict
    tuple val(patient), val(sample)

    output:
    tuple val(meta), path("*_filtered.vcf")      , emit: vcf

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${mutect_vcf.baseName}${options.suffix}" : "${mutect_vcf.baseName}"

    """
    gatk FilterMutectCalls \\
        --variant ${mutect_vcf} \\
        --stats ${mutect_vcf_stats} \\
        --output Mutect.vcf \\
        --reference ${fasta}

    mutect2_filter.py Mutect.vcf ${prefix}.vcf ${sample[1]} ${sample[0]}
    """
}
