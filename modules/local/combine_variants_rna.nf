// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process COMBINE_VARIANTS {
    tag "Combine_variants"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['patient']) }

    input:
    path(fasta)
    path(fai)
    path(dict)
    tuple val(meta), path(haplotypecaller_vcf_filtered)
    tuple val(meta), path(varscan_vcf)

    output:
    tuple val(meta), path("*_combined_calls.vcf"), path("*_combined_calls.vcf.idx")      , emit: vcf
    path "*.version.txt"                                                                 , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    gatk3 -T CombineVariants \\
        -R $fasta \\
        -V:varscan $varscan_snv_filtered \\
        -V:HaplotypeCaller $somaticsniper_filtered \\
        -o ${prefix}.vcf \\
        $options.args \\
        --num_threads $task.cpus

    sed -i 's/${meta.id}.HaplotypeCaller/HaplotypeCaller/g' ${prefix}.vcf

    sed -i 's/Sample1.varscan/varscan/g' combined_calls.vcf

    echo \$(gatk3 -T CombineVariants --version) > ${software}.version.txt
    """
}
