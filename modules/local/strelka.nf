// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process STRELKA {
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
    tuple val(meta), path("*.somatic_indels.vcf.gz")    , emit: vcf_indels
    tuple val(meta), path("*.somatic_indels.vcf.gz.tbi"), emit: vcf_indels_tbi
    tuple val(meta), path("*.somatic_snvs.vcf.gz")      , emit: vcf_snvs
    tuple val(meta), path("*.somatic_snvs.vcf.gz.tbi")  , emit: vcf_snvs_tbi
    path "*.version.txt"                                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def options_target_bed = useBed ? "--exome --callRegions call_targets.bed.gz" : ""
    def beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${target_bed} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""

    """
    ${beforeScript}
    configureStrelkaSomaticWorkflow.py \\
        --tumor $bam_tumor \\
        --normal $bam_normal \\
        --referenceFasta $fasta \\
        $options_target_bed \\
        $options.args \\
        --runDir strelka

    python2 strelka/runWorkflow.py -m local -j $task.cpus
    mv strelka/results/variants/somatic.indels.vcf.gz     ${prefix}.somatic_indels.vcf.gz
    mv strelka/results/variants/somatic.indels.vcf.gz.tbi ${prefix}.somatic_indels.vcf.gz.tbi
    mv strelka/results/variants/somatic.snvs.vcf.gz       ${prefix}.somatic_snvs.vcf.gz
    mv strelka/results/variants/somatic.snvs.vcf.gz.tbi   ${prefix}.somatic_snvs.vcf.gz.tbi

    echo \$( configureStrelkaSomaticWorkflow.py --version ) > ${software}.version.txt
    """
}
