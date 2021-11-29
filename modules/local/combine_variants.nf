// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COMBINE_VARIANTS {
    tag "Combine_variants"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    path(fasta)
    path(fai)
    path(dict)
    tuple val(meta), path(mutect2_filtered)
    tuple val(meta), path(varscan_indel_filtered)
    tuple val(meta), path(varscan_snv_filtered)
    tuple val(meta), path(strelka_indel_filtered)
    tuple val(meta), path(strelka_snv_filtered)
    tuple val(meta), path(somaticsniper_filtered)

    output:
    tuple val(meta), path("*_combined_calls.vcf"), path("*_combined_calls.vcf.idx")      , emit: vcf

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    gatk3 -T CombineVariants \\
        -R $fasta \\
        -V:varscan_indel $varscan_indel_filtered \\
        -V:varscan $varscan_snv_filtered \\
        -V:strelka_indel $strelka_indel_filtered \\
        -V:strelka $strelka_snv_filtered \\
        -V:mutect $mutect2_filtered \\
        -V:somaticsniper $somaticsniper_filtered \\
        -o ${prefix}.vcf \\
        $options.args \\
        --num_threads $task.cpus
    """
}
