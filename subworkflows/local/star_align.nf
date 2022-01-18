//
// Alignment with STAR
//

params.star_options          = [:]
params.samtools_index_options = [:]

include { STAR_ALIGN }     from '../../modules/local/star_align'       addParams( options: params.star_options )
include { SAMTOOLS_INDEX } from '../../modules/local/samtools_index'   addParams( options: params.samtools_index_options )

workflow ALIGN_STAR {
    take:
    reads // channel: [ val(meta), [ reads ] ]
    index // channel: /path/to/star/index/
    gtf   // channel: /path/to/genome.gtf

    main:

    ch_versions = Channel.empty()

    //
    // Map reads with STAR
    //
    STAR_ALIGN ( reads, index, gtf )
    ch_versions = ch_versions.mix(STAR_ALIGN.out.versions.first())

    bam = STAR_ALIGN.out.bam

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    SAMTOOLS_INDEX ( STAR_ALIGN.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    bam_indexed = bam.join(SAMTOOLS_INDEX.out.bai)

    emit:
    bam            = bam_indexed            // channel: [ val(meta), bam]
    versions       = ch_versions            // channel: [ versions.yml ]
}
