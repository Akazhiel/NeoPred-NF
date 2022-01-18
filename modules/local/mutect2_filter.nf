// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MUTECT2_FILTER {
    tag "$meta.id"
    label 'process_low'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['patient']) }

    input:
    tuple val(meta), path(mutect_vcf), path(mutect_vcf_stats)
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("${prefix}.vcf"), emit: vcf

    script:
    def software = getSoftwareName(task.process)
    prefix   = options.suffix ? "${mutect_vcf.baseName}${options.suffix}" : "${mutect_vcf.baseName}"

    """
    gatk FilterMutectCalls \\
        --variant ${mutect_vcf} \\
        --stats ${mutect_vcf_stats} \\
        --output Mutect2_${meta.id}.vcf \\
        --reference ${fasta}

    mutect2_filter.py Mutect2_${meta.id}.vcf ${prefix}.vcf ${meta.tumor} ${meta.normal}
    """
}
