//
// SOMATIC VARIANT CALLING
//

params.varscan_options                = [:]
params.haplotypecaller_options        = [:]

include { VARSCAN }                                      from '../../modules/local/varscan_rna'                  addParams(options: params.varscan_options)
include { HAPLOTYPECALLER }                              from '../../modules/local/gatk4/haplotypecaller'        addParams(options: params.haplotypecaller_options)

workflow RNA_VARIANT_CALLING {
    take:
        bam                   // channel: [mandatory] bam
        dbsnp                 // channel: [mandatory] dbsnp
        dbsnp_tbi             // channel: [mandatory] dbsnp_tbi
        dict                  // channel: [mandatory] dict
        fai                   // channel: [mandatory] fai
        fasta                 // channel: [mandatory] fasta

    main:

    haplotypecaller_vcf     = Channel.empty()
    varscan_vcf             = Channel.empty()
    tool_versions           = Channel.empty()


    HAPLOTYPECALLER (
        bam,
        fasta,
        fai,
        dict,
        dbsnp,
        dbsnp_tbi
    )

    haplotypecaller_vcf = HAPLOTYPECALLER.out.vcf
    tool_versions       = tool_versions.mix(HAPLOTYPECALLER.out.version)

    VARSCAN (
        bam,
        fasta,
        fai
    )

    varscan_vcf   = VARSCAN.out.vcf
    tool_versions = tool_versions.mix(VARSCAN.out.version)

    emit:
        haplotypecaller_vcf    = haplotypecaller_vcf
        varscan_vcf            = varscan_vcf
        version                = tool_versions
}
