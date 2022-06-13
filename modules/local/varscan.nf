// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process VARSCAN {
    tag "$meta.id"
    label 'process_high'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    input:
    tuple val(meta), path(bam_normal), path(bai_normal), path(bam_tumor), path(bai_tumor)
    path  fasta
    path  fai
    path  target_bed
    val   useBed

    output:
    tuple val(meta), path("*.snp.vcf")  , emit: vcf_snvs
    tuple val(meta), path("*.indel.vcf"), emit: vcf_indels
    path "*.version.txt"                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def options_target_bed = useBed ? "--positions ${target_bed}" : ""

    """

    samtools mpileup $options.args $options_target_bed -f ${fasta} ${bam_normal} > Normal.pileup
    samtools mpileup $options.args $options_target_bed -f ${fasta} ${bam_tumor} > Tumor.pileup

    varscan somatic    \\
        Normal.pileup  \\
        Tumor.pileup   \\
        ${prefix}      \\
        $options.args2
        
    echo \$(varscan 2>&1) | sed -e 's/^.*VarScan v//g; s/ \\*.*\$//' > ${software}.version.txt
    """
}
