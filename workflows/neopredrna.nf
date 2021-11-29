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
                            params.dbsnp,
                            params.dbsnp_index,
                            params.known_indels,
                            params.known_indels_index,
                            params.hla_reference_rna,
                            params.gtf
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
star_index              = Channel.fromPath(params.star_index).collect()
known_indels            = Channel.fromPath(params.known_indels).collect()
known_indels_index      = Channel.fromPath(params.known_indels_index).collect()
dbsnp                   = Channel.fromPath(params.dbsnp).collect()
dbsnp_tbi               = Channel.fromPath(params.dbsnp_index).collect()
hla_rna_fasta           = Channel.fromPath(params.hla_reference_rna).collect()
dict                    = Channel.fromPath(params.dict).collect()
fasta                   = Channel.fromPath(params.fasta).collect()
fasta_fai               = Channel.fromPath(params.fasta_fai).collect()
gtf                     = Channel.fromPath(params.gtf).collect()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

include { TRIMGALORE }              from '../modules/local/trimgalore'                      addParams( options: modules['trimgalore'] )

include { STAR_ALIGN }              from '../modules/local/star_align'                      addParams( options: modules['star_align'] )

include { MARKDUPLICATES }          from '../subworkflows/local/markduplicates'             addParams(
    markduplicates_options:             modules['markduplicates'],
    markduplicatesspark_options:        modules['markduplicatesspark']
)

include { PREPARE_RECALIBRATION } from '../subworkflows/local/prepare_recalibration'        addParams(
    baserecalibrator_options:          modules['baserecalibrator'],
    baserecalibrator_spark_options:    modules['baserecalibrator_spark']
)

include { RECALIBRATE } from '../subworkflows/local/recalibrate'                            addParams(
    applybqsr_options:                 modules['applybqsr'],
    applybqsr_spark_options:           modules['applybqsr_spark'],
    qualimap_bamqc_options:            modules['qualimap_bamqc_recalibrate']
)

include { HLATYPING } from '../subworkflows/local/hlatyping'                                addParams(
    yara_index:                        modules['yara_index'],
    yara_mapping:                      modules['yara_map'],
    optitype:                          modules['optitype']
)

include { RNA_VARIANT_CALLING } from '../subworkflows/local/rna_variant_calling'            addParams(

)

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
// def multiqc_report = []

workflow NEOPRED_RNA {

    ch_software_versions = Channel.empty()

    // TRIM READS WITH CUTADAPT

    TRIMGALORE (
        input_sample
    )

    trimmed_reads = TRIMGALORE.out.reads

    ch_software_versions = ch_software_versions.mix(TRIMGALORE.out.version.first().ifEmpty(null))

    // ALIGN READS WITH STAR

    STAR_ALIGN (
        trimmed_reads,
        star_index,
        gtf
    )

    bam = STAR_ALIGN.out.bam

    ch_software_versions = ch_software_versions.mix(STAR_ALIGN.out.version.ifEmpty(null))

    // MARK DUPLICATES

    MARKDUPLICATES (
        bam,
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

    table_bqsr = PREPARE_RECALIBRATION.out.table_bqsr
    bam_applybqsr = bam_markduplicates.join(table_bqsr)
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
        hla_rna_fasta
    )

    ch_software_versions = ch_software_versions.mix(HLATYPING.out.version.ifEmpty(null))

    //
    // VARIANT CALLING
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
        if (row.type == "rna") {
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