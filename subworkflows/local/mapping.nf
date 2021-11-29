//
// MAPPING
//

params.bwamem_options = [:]
params.bwamem_tumor_options   = [:]
params.samtools_index_options = [:]

include { BWA_MEM as BWAMEM_N }       from '../../modules/local/bwamem'                    addParams(options: params.bwamem_options)
include { BWA_MEM as BWAMEM_T }       from '../../modules/local/bwamem'                    addParams(options: params.bwamem_tumor_options)
include { SAMTOOLS_INDEX }            from '../../modules/local/samtools_index'            addParams(options: params.samtools_index_options)

workflow MAPPING {
    take:
        reads_input         // channel: [mandatory] meta, reads_input
        bwa                 // channel: [mandatory] bwa

    main:

    bam_indexed = Channel.empty()

    reads_input_split = reads_input

    // If meta.status is 1, then sample is tumor
    // else, (even is no meta.status exist) sample is normal
    reads_input_split.branch{ meta, reads ->
        tumor:  meta.status == 1
        normal: true
    }.set{reads_input_status}

    bam_bwamem       = Channel.empty()
    bam_bwamem_n     = Channel.empty()
    bam_bwamem_t     = Channel.empty()
    tool_versions    = Channel.empty()

    BWAMEM_N(reads_input_status.normal, bwa)
    BWAMEM_T(reads_input_status.tumor, bwa)

    bam_bwamem_n = BWAMEM_N.out.bam
    bam_bwamem_t = BWAMEM_T.out.bam
    bam_bwamem   = bam_bwamem_n.mix(bam_bwamem_t)

    bwamem_n_version = BWAMEM_N.out.version
    bwamem_t_version = BWAMEM_T.out.version

    bwamem_version = bwamem_n_version.mix(bwamem_t_version).first()

    tool_versions = tool_versions.mix(bwamem_version)


    SAMTOOLS_INDEX(bam_bwamem)
    bam_indexed = bam_bwamem.join(SAMTOOLS_INDEX.out.bai)

    emit:
        bam         = bam_bwamem
        bam_indexed = bam_indexed
        versions    = tool_versions
}
