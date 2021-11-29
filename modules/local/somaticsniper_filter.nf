// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SOMATICSNIPER_FILTER {
    tag "Somaticsniper_filtering"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(somaticsniper_vcf)

    output:
    tuple val(meta), path("*_filtered.vcf")      , emit: vcf

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${somaticsniper_vcf.baseName}${options.suffix}" : "${somaticsniper_vcf.baseName}"

    """
    somaticsniper_filter.py ${somaticsniper_vcf}

    awk '{if (\$1 ~ /#/) {print} else if (\$4 != \$5) {gsub(/W|K|B|Y|D|H|V|R|S|M/,"N",\$4); OFS="\t"; print}}' tmp_ss.vcf > ${prefix}.vcf
    """
}
