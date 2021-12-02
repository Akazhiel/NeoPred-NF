/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowNeoprednf.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
                            params.input,
                            params.fasta,
                            params.fasta_fai,
                            params.dict,
                            params.bwa,
                            params.germline,
                            params.germline_index,
                            params.pon,
                            params.pon_index,
                            params.dbsnp,
                            params.dbsnp_index,
                            params.known_indels,
                            params.known_indels_index,
                            params.hla_reference_dna
                         ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()

input_sample            = extract_csv(ch_input)
bwa_index               = Channel.fromPath(params.bwa).collect()
known_indels            = Channel.fromPath(params.known_indels).collect()
known_indels_index      = Channel.fromPath(params.known_indels_index).collect()
dbsnp                   = Channel.fromPath(params.dbsnp).collect()
dbsnp_tbi               = Channel.fromPath(params.dbsnp_index).collect()
target_bed              = params.target_bed        ? Channel.fromPath(params.target_bed).collect()        : ch_dummy_file
hla_dna_fasta           = Channel.fromPath(params.hla_reference_dna).collect()
panel_of_normals        = Channel.fromPath(params.pon).collect()
panel_of_normals_tbi    = Channel.fromPath(params.pon_index).collect()
dict                    = Channel.fromPath(params.dict).collect()
fasta                   = Channel.fromPath(params.fasta).collect()
fasta_fai               = Channel.fromPath(params.fasta_fai).collect()
germline_resource       = Channel.fromPath(params.germline).collect()
germline_resource_tbi   = Channel.fromPath(params.germline_index).collect()
vep_cache_version       = params.vep_cache_version ?: Channel.empty()
vep_cache               = params.vep_cache ? Channel.fromPath(params.vep_cache).collect()                 : ch_dummy_file
vep_genome              = params.vep_genome ?: Channel.empty()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

// def multiqc_options   = modules['multiqc']
// multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//

//include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

include { FASTQC  }                 from '../modules/local/fastqc'                          addParams( options: modules['fastqc'] )

include { TRIMGALORE }              from '../modules/local/trimgalore'                      addParams( options: modules['trimgalore'] )

include { MAPPING } from '../subworkflows/local/mapping' addParams(
    bwamem_options:                     modules['bwamem'],
    bwamem_tumor_options:               modules['bwamem_tumor'],
    samtools_index_options:             modules['samtools_index_bam']
)

include { MARKDUPLICATES }          from '../subworkflows/local/markduplicates'             addParams(
    markduplicates_options: modules['markduplicates'],
    markduplicatesspark_options: modules['markduplicatesspark']
)

include { PREPARE_RECALIBRATION } from '../subworkflows/local/prepare_recalibration' addParams(
    baserecalibrator_options:          modules['baserecalibrator'],
    baserecalibrator_spark_options:    modules['baserecalibrator_spark']
)

include { RECALIBRATE } from '../subworkflows/local/recalibrate' addParams(
    applybqsr_options:                 modules['applybqsr'],
    applybqsr_spark_options:           modules['applybqsr_spark'],
    qualimap_bamqc_options:            modules['qualimap_bamqc_recalibrate']
)

include { HLATYPING } from '../subworkflows/local/hlatyping' addParams(
    yara_index:                        modules['yara_index'],
    yara_mapping:                      modules['yara_map'],
    optitype:                          modules['optitype']
)

include { PAIR_VARIANT_CALLING } from '../subworkflows/local/somatic_variant_calling' addParams(
    mutect2_somatic_options:            modules['mutect2_somatic'],
    strelka_options:                    modules['strelka_somatic'],
    somaticsniper_options:              modules['somaticsniper'],
    varscan_options:                    modules['varscan']
)

include { FILTER_VARIANTS } from '../subworkflows/local/filter_variants' addParams(
    strelka_filter:                     modules['strelka_filter'],
    mutect2_filter:                     modules['mutect2_filter'],
    varscan_filter:                     modules['varscan_filter'],
    somaticsniper_filter:               modules['somaticsniper_filter']
)

include { COMBINE_VARIANTS } from '../modules/local/combine_variants' addParams(options:  modules['combine_variants'])

include { VEP } from '../modules/local/vep_annotate' addParams(options: modules['vep'])

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
// def multiqc_report = []

workflow NEOPRED_DNA {

    ch_software_versions = Channel.empty()

    known_sites     = dbsnp.concat(known_indels).collect()
    known_sites_tbi = dbsnp_tbi.concat(known_indels_index).collect()

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        input_sample
    )

    ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    //
    // MODULE: Run Trimgalore
    //

    TRIMGALORE (
        input_sample
    )

    trimmed_reads = TRIMGALORE.out.reads

    ch_software_versions = ch_software_versions.mix(TRIMGALORE.out.version.ifEmpty(null))

    MAPPING (
        trimmed_reads,
        bwa_index
    )

    bam_mapped  = MAPPING.out.bam
    bam_indexed = MAPPING.out.bam_indexed

    ch_software_versions = ch_software_versions.mix(MAPPING.out.versions.ifEmpty(null))

    // trimmed_reads.branch{ meta, reads ->
    //     tumor:  meta.status == 1
    //     normal: true
    // }.set{reads_input_status}

    // BWA_MEM_N (
    //     reads_input_status.normal,
    //     bwa_index
    // )

    // BWA_MEM_T (
    //     reads_input_status.tumor,
    //     bwa_index
    // )

    // bam_bwamem_n = BWA_MEM_N.out.bam
    // bam_bwamem_t = BWA_MEM_T.out.bam
    // bam_bwamem = bam_bwamem_n.mix(bam_bwamem_t)

    // ch_software_versions = ch_software_versions.mix(BWA_MEM_N.out.version.ifEmpty(null))

    //
    // Mark duplicates
    //

    MARKDUPLICATES (
        bam_mapped,
        params.use_metrics,
        params.use_gatk_spark,
        fasta,
        fasta_fai,
        dict
    )

    bam_markduplicates = MARKDUPLICATES.out.bam

    ch_software_versions = ch_software_versions.mix(MARKDUPLICATES.out.versions.ifEmpty(null))

    //
    // Prepare Recalibration
    //

    PREPARE_RECALIBRATION (
        bam_markduplicates,
        params.use_gatk_spark,
        dict,
        fasta_fai,
        fasta,
        known_sites,
        known_sites_tbi,
    )

    table_bqsr           = PREPARE_RECALIBRATION.out.table_bqsr
    bam_applybqsr        = bam_markduplicates.join(table_bqsr)
    ch_software_versions = ch_software_versions.mix(PREPARE_RECALIBRATION.out.version.ifEmpty(null))

    //
    // Apply Recalibration
    //

    RECALIBRATE (
        params.use_gatk_spark,
        bam_applybqsr,
        dict,
        fasta_fai,
        fasta,
        target_bed
    )

    bam_bqsr = RECALIBRATE.out.bam
    ch_software_versions = ch_software_versions.mix(RECALIBRATE.out.version.ifEmpty(null))

    //
    // HLATYPING
    //

    HLATYPING (
        bam_bqsr,
        hla_dna_fasta
    )

    hla = HLATYPING.out.hla
    ch_software_versions = ch_software_versions.mix(HLATYPING.out.version.ifEmpty(null))

    //
    // SOMATIC VARIANT CALLING
    //

    PAIR_VARIANT_CALLING (
        bam_bqsr,
        dbsnp,
        dbsnp_tbi,
        dict,
        fasta_fai,
        fasta,
        target_bed,
        germline_resource,
        germline_resource_tbi,
        panel_of_normals,
        panel_of_normals_tbi
    )

    mutect2_vcf       = PAIR_VARIANT_CALLING.out.mutect2_vcf
    mutect2_vcf_stats = PAIR_VARIANT_CALLING.out.mutect2_vcf_stats
    varscan_vcf       = PAIR_VARIANT_CALLING.out.varscan_vcf
    somaticsniper_vcf = PAIR_VARIANT_CALLING.out.somaticsniper_vcf
    strelka_vcf       = PAIR_VARIANT_CALLING.out.strelka_vcf

    ch_software_versions = ch_software_versions.mix(PAIR_VARIANT_CALLING.out.version.ifEmpty(null))

    //
    // SOMATIC FILTERING
    //

    FILTER_VARIANTS (
        bam_bqsr,
        fasta,
        fasta_fai,
        dict,
        varscan_vcf,
        somaticsniper_vcf,
        strelka_vcf,
        mutect2_vcf,
        mutect2_vcf_stats
    )

    mutect2_filtered        = FILTER_VARIANTS.out.mutect2_vcf
    varscan_indel_filtered  = FILTER_VARIANTS.out.varscan_indel_vcf
    varscan_snv_filtered    = FILTER_VARIANTS.out.varscan_snv_vcf
    strelka_indel_filtered  = FILTER_VARIANTS.out.strelka_indel_vcf
    strelka_snv_filtered    = FILTER_VARIANTS.out.strelka_snv_vcf
    somaticsniper_filtered  = FILTER_VARIANTS.out.somaticsniper_vcf

    //
    // MERGE VARIANTS
    //

    COMBINE_VARIANTS (
        fasta,
        fasta_fai,
        dict,
        mutect2_filtered,
        varscan_indel_filtered,
        varscan_snv_filtered,
        strelka_indel_filtered,
        strelka_snv_filtered,
        somaticsniper_filtered
    )

    merged_vcf = COMBINE_VARIANTS.out.vcf
    merged_vcf.dump()
    ch_software_versions = ch_software_versions.mix(COMBINE_VARIANTS.out.version.ifEmpty(null))


    //
    //  ANNOTATE VARIANTS
    //

    VEP (
        merged_vcf,
        fasta,
        vep_cache,
        vep_cache_version,
        vep_genome
    )

    annotated_vcf = VEP.out.vcf.collect().dump()
    ch_software_versions = ch_software_versions.mix(VEP.out.version.ifEmpty(null))

    //
    // MODULE: Pipeline reporting
    //

    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // MODULE: MultiQC
    //
    // workflow_summary    = WorkflowNeoprednf.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    // MULTIQC (
    //     ch_multiqc_files.collect()
    // )
    // multiqc_report       = MULTIQC.out.report.toList()
    // ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))

    emit:
    vep_vcf = annotated_vcf
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

def extract_csv(csv_file) {
    Channel.from(csv_file).splitCsv(header: true)
        //Retrieves number of lanes by grouping together by patient and sample and counting how many entries there are for this combination
        .map{ row ->
            if (!(row.patient && row.sample)) log.warn "Missing or unknown field in csv file header"
            [[row.patient.toString(), row.sample.toString()], row]
        }.groupTuple()
        .map{ meta, rows ->
            size = rows.size()
            [rows, size]
        }.transpose()
        .map{ row, numLanes -> //from here do the usual thing for csv parsing
        def meta = [:]

        //TODO since it is mandatory: error/warning if not present?
        // Meta data to identify samplesheet
        // Both patient and sample are mandatory
        // Several sample can belong to the same patient
        // Sample should be unique for the patient
        if (row.type == "dna") {
            if (row.patient) meta.patient = row.patient.toString()
            if (row.sample)  meta.sample  = row.sample.toString()

            // If no gender specified, gender is not considered
            // gender is only mandatory for somatic CNV
            if (row.gender) meta.gender = row.gender.toString()
            else meta.gender = "NA"

            // If no status specified, sample is assumed normal
            if (row.status) meta.status = row.status.toInteger()
            else meta.status = 0

            // Add if sample type is DNA or RNA
            if (row.type) meta.seqtype = row.type.toString()
            else meta.seqtype = "NA"

            // mapping with fastq
            if (row.fastq_2) {
                meta.id         = "${row.patient}_${row.sample}".toString()
                def fastq_1     = file(row.fastq_1, checkIfExists: true)
                def fastq_2     = file(row.fastq_2, checkIfExists: true)
                def CN          = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ''
                def read_group  = "\"@RG\\tID:${meta.id}\\t${CN}PU:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}\\tPL:ILLUMINA\""
                meta.numLanes = numLanes.toInteger()
                meta.read_group = read_group.toString()
                return [meta, [fastq_1, fastq_2]]
            } else {
                log.warn "Missing or unknown field in csv file header"
            }
        }
    }
}

/*
========================================================================================
    THE END
========================================================================================
*/
