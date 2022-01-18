//
// PREPARE RECALIBRATION
//

params.baserecalibrator_options        = [:]
params.baserecalibrator_spark_options  = [:]

include { GATK4_BASERECALIBRATOR  as BASERECALIBRATOR             } from '../../modules/local/gatk4/baserecalibrator'      addParams(options: params.baserecalibrator_options)
include { GATK4_BASERECALIBRATOR_SPARK  as BASERECALIBRATOR_SPARK } from '../../modules/local/gatk4/baserecalibratorspark' addParams(options: params.baserecalibrator_spark_options)

workflow PREPARE_RECALIBRATION {
    take:
        bam  // channel: [mandatory] bam_markduplicates
        use_gatk_spark      //   value: [mandatory] use gatk spark
        dict                // channel: [mandatory] dict
        fai                 // channel: [mandatory] fai
        fasta               // channel: [mandatory] fasta
        known_sites         // channel: [optional]  known_sites
        known_sites_tbi     // channel: [optional]  known_sites_tbi
        target_bed

    main:

    tool_versions         = Channel.empty()

    if (use_gatk_spark) {
        BASERECALIBRATOR_SPARK(bam, fasta, fai, dict, known_sites, known_sites_tbi)
        table_baserecalibrator = BASERECALIBRATOR_SPARK.out.table
        tool_versions          = BASERECALIBRATOR_SPARK.out.version
    } else {
        BASERECALIBRATOR(bam, fasta, fai, dict, known_sites, known_sites_tbi, target_bed)
        table_baserecalibrator = BASERECALIBRATOR.out.table
        tool_versions          = BASERECALIBRATOR.out.version
    }

    //STEP 3.5: MERGING RECALIBRATION TABLES
    table_baserecalibrator.map { meta, table ->
        [meta, table]
    }.set{table_bqsr}

    emit:
        table_bqsr = table_bqsr
        version    = tool_versions
}
