#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GATK4_BASERECALIBRATOR as BASERECALIBRATOR } from '../modules/local/gatk4/baserecalibrator' addParams(options: [publish_dir: 'preprocessing', publish_by_meta: true, publish_files: ['recal.table': 'recal_table'] ])

workflow test_gatk_recalibrator {

    input = [
        [patient:'GX001', sample:'FrTu', gender:'NA', status:1, id:'GX001_FrTu', numLanes:1, read_group:'"@RG\tID:FrTu\tPU:FrTu\tSM:FrTu\tLB:FrTu\tPL:ILLUMINA"'], file('/home/jonathan/Nano_docker/HLA_NF/nf-core-neoprednf/work/8e/e76b1f952851f1f7b6d9b9baa8aff5//GX001_FrTu_md.bam', checkIfExists: true),
        file('/home/jonathan/Nano_docker/HLA_NF/nf-core-neoprednf/work/8e/e76b1f952851f1f7b6d9b9baa8aff5//GX001_FrTu_md.bam.bai', checkIfExists: true)
    ]

    fasta = file('/media/AGROS/References/hg19/Homo_sapiens_assembly19.fasta', checkIfExists: true)
    fai = file('/media/AGROS/References/hg19/Homo_sapiens_assembly19.fasta.fai', checkIfExists: true)
    dict = file('/media/AGROS/References/hg19/Homo_sapiens_assembly19.dict', checkIfExists: true)
    sites = file('/media/AGROS/References/hg19/dbsnp_138.b37.vcf.gz', checkIfExists: true)
    sites_tbi = file('/media/AGROS/References/hg19/dbsnp_138.b37.vcf.gz.tbi', checkIfExists: true)

    println "${input}"

    BASERECALIBRATOR(
        input,
        fasta,
        fai,
        dict,
        sites,
        sites_tbi
    )
}

workflow {
    test_gatk_recalibrator ()
}
