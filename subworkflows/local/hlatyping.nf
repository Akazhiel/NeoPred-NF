//
// MAPPING
//

params.yara_index   = [:]
params.yara_mapping = [:]
params.optitype     = [:]

include { YARA_INDEX }       from '../../modules/local/yaraindex'           addParams(options: params.yara_index)
include { YARA_MAPPER }      from '../../modules/local/yaramap'             addParams(options: params.yara_mapping)
include { OPTITYPE }         from '../../modules/local/optitype'            addParams(options: params.optitype)

workflow HLATYPING {
    take:
        bam                       // channel: [mandatory] meta, bam recalibrated
        hla_fasta                 // channel: [mandatory] bwa
        input

    main:

    tool_versions = Channel.empty()

    // Index hla reference fasta

    YARA_INDEX(hla_fasta, input)
    hla_index = YARA_INDEX.out.index
    tool_versions = tool_versions.mix(YARA_INDEX.out.version)

    // Map the reads to the hla reference

    YARA_MAPPER(bam, hla_index)
    hla_bam = YARA_MAPPER.out.bam
    tool_versions = tool_versions.mix(YARA_MAPPER.out.version)
    // Perform hlatyping with OptiType

    OPTITYPE(hla_bam)
    tool_versions = tool_versions.mix(OPTITYPE.out.version)
    hlatype = OPTITYPE.out.output

    emit:
        hla         = hlatype
        version     = tool_versions
}
