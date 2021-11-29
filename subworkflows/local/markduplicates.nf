//
// MARKDUPLICATES AND/OR QC after mapping
//

params.markduplicates_options            = [:]
params.markduplicatesspark_options       = [:]

include { GATK4_MARKDUPLICATES }                          from '../../modules/local/gatk4/markduplicates'             addParams(options: params.markduplicates_options)
include { GATK4_MARKDUPLICATES_SPARK }                    from '../../modules/local/gatk4/markduplicatesspark'        addParams(options: params.markduplicatesspark_options)

workflow MARKDUPLICATES {
    take:
        bam_mapped          // channel: [mandatory] meta, bam
        save_metrics
        use_gatk_spark      // value: [mandatory] use gatk spark
        fasta               // channel: [mandatory] fasta
        fai                 // channel: [mandatory] fai
        dict                // channel: [mandatory] dict

    main:

    report_markduplicates = Channel.empty()
    tool_versions         = Channel.empty()

    if (use_gatk_spark) {
        //If BAMQC should be run on MD output, then don't use MDSpark to convert to cram, but use bam output instead
        GATK4_MARKDUPLICATES_SPARK(bam_mapped, fasta, fai, dict)
        bam_markduplicates = GATK4_MARKDUPLICATES_SPARK.out.output
        tool_versions = tool_versions.mix(GATK4_MARKDUPLICATES_SPARK.out.version)
    } else {
        GATK4_MARKDUPLICATES(bam_mapped, save_metrics)
        report_markduplicates = GATK4_MARKDUPLICATES.out.metrics
        bam_markduplicates    = GATK4_MARKDUPLICATES.out.bam
        tool_versions = tool_versions.mix(GATK4_MARKDUPLICATES.out.version)
    }

    emit:
        bam         = bam_markduplicates
        versions    = tool_versions
        // qc      = qc_reports
}
