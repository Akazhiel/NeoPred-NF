//
// SOMATIC FILTERING
//

params.strelka_filter       = [:]
params.mutect2_filter       = [:]
params.varscan_filter       = [:]
params.somaticsniper_filter = [:]

include { STRELKA_FILTER_INDEL }        from '../../modules/local/strelka_filter_indel'         addParams(options: params.strelka_filter)
include { STRELKA_FILTER_SNV }          from '../../modules/local/strelka_filter_snv'           addParams(options: params.strelka_filter)
include { SOMATICSNIPER_FILTER }        from '../../modules/local/somaticsniper_filter'         addParams(options: params.somaticsniper_filter)
include { MUTECT2_FILTER }              from '../../modules/local/mutect2_filter'               addParams(options: params.mutect2_filter)
include { VARSCAN_FILTER }              from '../../modules/local/varscan_filter'               addParams(options: params.varscan_filter)

workflow FILTER_VARIANTS {
    take:
        bam
        fasta
        fasta_fai
        dict
        varscan_vcf
        somaticsniper_vcf
        strelka_vcf
        mutect2_vcf
        mutect2_vcf_stats

    main:

    tool_versions = Channel.empty()

    bam.map{ meta, bam, bai ->
        patient = meta.patient
        sample  = meta.id
        status  = meta.status
        [patient, sample]
    }.groupTuple(sort: true).set{ bam }

    strelka_vcf.branch{ meta, vcfs ->
        indel:  vcfs =~ /\\*_indels.vcf.gz/
        snv:    vcfs =~ /\\*_snvs.vcf.gz/
    }.set { strelka_vcf_split }

    STRELKA_FILTER_INDEL (
        strelka_vcf_split.indel
    )

    strelka_indel_vcf = STRELKA_FILTER_INDEL.out.vcf

    STRELKA_FILTER_SNV (
        strelka_vcf_split.snv
    )

    strelka_snv_vcf = STRELKA_FILTER_SNV.out.vcf

    SOMATICSNIPER_FILTER (
        somaticsniper_vcf
    )

    somaticsniper_vcf = SOMATICSNIPER_FILTER.out.vcf

    MUTECT2_FILTER (
        mutect2_vcf,
        mutect2_vcf_stats,
        fasta,
        fasta_fai,
        dict,
        bam
    )

    mutect2_vcf = MUTECT2_FILTER.out.vcf

    VARSCAN_FILTER (
        varscan_vcf
    )

    varscan_snv_vcf   = VARSCAN_FILTER.out.snv_vcf
    varscan_indel_vcf = VARSCAN_FILTER.out.indel_vcf

    emit:
        strelka_indel_vcf = strelka_indel_vcf
        strelka_snv_vcf   = strelka_snv_vcf
        somaticsniper_vcf = somaticsniper_vcf
        mutect2_vcf       = mutect2_vcf
        varscan_snv_vcf   = varscan_snv_vcf
        varscan_indel_vcf = varscan_indel_vcf
}
