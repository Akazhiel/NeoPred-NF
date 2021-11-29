//
// SOMATIC VARIANT CALLING
//

params.varscan_options                = [:]
params.strelka_options                = [:]
params.mutect2_somatic_options        = [:]
params.somaticsniper_options          = [:]

include { VARSCAN }                                      from '../../modules/local/varscan'                  addParams(options: params.varscan_options)
include { HAPLOTYPECALLERCDS790¡ }                                from '../../modules/local/somaticsniper'            addParams(options: params.somaticsniper_options)

workflow PAIR_VARIANT_CALLING {
    take:
        bam                   // channel: [mandatory] bam
        dbsnp                 // channel: [mandatory] dbsnp
        dbsnp_tbi             // channel: [mandatory] dbsnp_tbi
        dict                  // channel: [mandatory] dict
        fai                   // channel: [mandatory] fai
        fasta                 // channel: [mandatory] fasta
        target_bed            // channel: [optional]  target_bed
        germline_resource     // channel: [optional]  germline_resource
        germline_resource_tbi // channel: [optional]  germline_resource_tbi
        panel_of_normals      // channel: [optional]  panel_of_normals
        panel_of_normals_tbi  // channel: [optional]  panel_of_normals_tbi

    main:

    bam.map{ meta, bam, bai ->
        patient = meta.patient
        sample  = meta.id
        gender  = meta.gender
        status  = meta.status
        [patient, sample, gender, status, bam, bai]
    }.branch{
        normal: it[3] == 0
        tumor:  it[3] == 1
    }.set{ bam_to_cross }

    bam_pair = bam_to_cross.normal.cross(bam_to_cross.tumor).map { normal, tumor ->
        def meta = [:]
        meta.patient = normal[0]
        meta.normal  = normal[1]
        meta.tumor   = tumor[1]
        meta.gender  = normal[2]
        meta.id      = "${meta.tumor}_vs_${meta.normal}".toString()

        [meta, normal[4], normal[5], tumor[4], tumor[5]]
    }

    haplotypecaller_vcf     = Channel.empty()
    varscan_vcf             = Channel.empty()
    tool_versions           = Channel.empty()


    SOMATICSNIPER (
        bam_pair,
        fasta,
        fai
    )

    somaticsniper_vcf = SOMATICSNIPER.out.vcf
    tool_versions = tool_versions.mix(SOMATICSNIPER.out.version)

    VARSCAN (
        bam_pair,
        fasta,
        fai,
        target_bed,
        params.target_bed
    )

    varscan_snvs_vcf   = VARSCAN.out.vcf

    varscan_vcf = varscan_vcf.mix(varscan_snvs_vcf, varscan_indels_vcf)
    tool_versions = tool_versions.mix(VARSCAN.out.version)


    MUTECT2 (
        bam_pair,
        panel_of_normals,
        panel_of_normals_tbi,
        dict,
        fasta,
        fai,
        germline_resource,
        germline_resource_tbi,
        target_bed,
        params.target_bed
    )

    mutect2_vcf = MUTECT2.out.vcf.collect()
    tool_versions = tool_versions.mix(MUTECT2.out.version)

    emit:
        mutect2_vcf            = mutect2_vcf
        mutect2_vcf_stats      = mutect2_vcf_stats
        strelka_vcf            = strelka_vcf
        varscan_vcf            = varscan_vcf
        somaticsniper_vcf      = somaticsniper_vcf
        version                = tool_versions
}
