// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process VARSCAN_FILTER {
    tag "Varscan_filtering"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['sample']) }

    input:
    tuple val(meta), path(varscan_vcf)

    output:
    tuple val(meta), path("*snp_filtered.vcf"), optional: true        , emit: snv_vcf
    tuple val(meta), path("*indel_filtered.vcf"), optional: true      , emit: indel_vcf

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${varscan_vcf.baseName}${options.suffix}" : "${varscan_vcf.baseName}"

    """
    varscan_filter.py ${varscan_vcf}

    awk '{if (\$1 ~ /#/) {print} else if (\$4 != \$5) {gsub(/W|K|B|Y|D|H|V|R|S|M/,"N",\$4); OFS="\t"; print}}' tmp_varscan.vcf > ${prefix}.vcf

    """
}
