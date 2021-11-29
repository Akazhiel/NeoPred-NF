//
// RECALIBRATE
//

params.applybqsr_options      = [:]
params.applybqsr_spark_options = [:]
params.qualimap_bamqc_options = [:]

include { GATK4_APPLYBQSR as APPLYBQSR }             from '../../modules/local/gatk4/applybqsr'      addParams(options: params.applybqsr_options)
include { GATK4_APPLYBQSR_SPARK as APPLYBQSR_SPARK } from '../../modules/local/gatk4/applybqsrspark' addParams(options: params.applybqsr_spark_options)
include { QUALIMAP_BAMQC }                           from '../../modules/local/qualimap'             addParams(options: params.qualimap_bamqc_options)

workflow RECALIBRATE {
    take:
        use_gatk_spark //   value: [mandatory] use gatk spark
        bam            // channel: [mandatory] bam
        dict           // channel: [mandatory] dict
        fai            // channel: [mandatory] fai
        fasta          // channel: [mandatory] fasta
        target_bed     // channel: [optional]  target_bed

    main:

    bam_recalibrated         = Channel.empty()
    tool_versions            = Channel.empty()

    if(use_gatk_spark){
        APPLYBQSR_SPARK(bam, fasta, fai, dict)
        bam_applybqsr = APPLYBQSR_SPARK.out.bam
        tool_versions = tool_versions.mix(APPLYBQSR_SPARK.out.version)
    }else{
        APPLYBQSR(bam, fasta, fai, dict)
        bam_applybqsr = APPLYBQSR.out.bam
        tool_versions = tool_versions.mix(APPLYBQSR.out.version)
    }

    QUALIMAP_BAMQC(bam_applybqsr, target_bed, params.target_bed)
    qualimap_bamqc = QUALIMAP_BAMQC.out.results
    tool_versions = tool_versions.mix(QUALIMAP_BAMQC.out.version)

    emit:
        bam     = bam_applybqsr
        qc      = qualimap_bamqc
        version = tool_versions
}
