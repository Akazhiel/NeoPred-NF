// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process MUTECT2 {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor)
    path pon
    path ponIndex
    path dict
    path fasta
    path fai
    path(germline_resource)
    path(germline_resource_tbi)
    path target_bed
    val useBed

    output:
    tuple val(meta), path("*.vcf"),       emit: vcf
    tuple val(meta), path("*.vcf.stats"), emit: vcf_stats
    path "*.version.txt"                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def intervalsOptions = useBed ? "-L ${target_bed}" : ""
    def softClippedOption = params.ignore_soft_clipped_bases ? "--dont-use-soft-clipped-bases true" : ""
    def PON = "--panel-of-normals ${pon}"

    """
    # Get raw calls
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        Mutect2 \
        -R ${fasta}\
        -I ${bam_tumor} -tumor ${meta.tumor} \
        -I ${bam_normal} -normal ${meta.normal} \
        ${intervalsOptions} \
        ${softClippedOption} \
        --germline-resource ${germline_resource} \
        ${PON} \
        -O ${prefix}.vcf
    echo \$(gatk Mutect2 --version 2>&1) | sed 's/^.*(GATK) v//; s/ HTSJDK.*\$//' > ${software}.version.txt
    """
}
