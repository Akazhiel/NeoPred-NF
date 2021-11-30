// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process VARSCAN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam), path(bai)
    path(fasta)
    path(fai)
    path(dbsnp)
    path(dbsnp_index)

    output:
    tuple val(meta), path("*.vcf")  , emit: vcf
    path "*.version.txt"            , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        HaplotypeCaller \
        --reference ${fasta}\
        --input ${bam} \
        --dbsnp ${dbsnp} \
        $options.args \
        --output ${prefix}.vcf

    echo \$(gatk HaplotypeCaller --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}
