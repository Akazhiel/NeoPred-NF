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
target_bed              = params.target_bed        ? Channel.fromPath(params.target_bed).collect()        : ch_dummy_file
dbsnp_tbi               = Channel.fromPath(params.dbsnp_index).collect()
hla_rna_fasta           = Channel.fromPath(params.hla_reference_rna).collect()
dict                    = Channel.fromPath(params.dict).collect()
fasta                   = Channel.fromPath(params.fasta).collect()
fasta_fai               = Channel.fromPath(params.fasta_fai).collect()
gtf                     = Channel.fromPath(params.gtf).collect()
vep_cache_version       = params.vep_cache_version ?: Channel.empty()
vep_cache               = params.vep_cache ? Channel.fromPath(params.vep_cache).collect()                 : []
vep_genome              = params.vep_genome ?: Channel.empty()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'RNA']] )

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

include { SPLITNCIGAR }             from '../modules/local/gatk4/splitreads'                addParams( options: modules['splitncigarreads'])

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
    varscan_options:                   modules['varscan_rna'],
    haplotypecaller_options:           modules['haplotypecaller']
)

include { HAPLOTYPECALLER_FILTER } from '../modules/local/gatk4/variantfiltration'          addParams( options: modules['haplotypecaller_filter'])

include { FEATURECOUNTS } from '../modules/local/featurecounts'                             addParams( options: modules['featurecounts'])

include { COMBINE_VARIANTS } from '../modules/local/combine_variants_rna'                   addParams( options: modules['combine_variants'])

include { VEP } from '../modules/local/vep_annotate'                                        addParams( options: modules['vep'])

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
// def multiqc_report = []

workflow NEOPRED_RNA {

    ch_software_versions = Channel.empty()

    known_sites     = dbsnp.concat(known_indels).collect()
    known_sites_tbi = dbsnp_tbi.concat(known_indels_index).collect()

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

    ch_software_versions = ch_software_versions.mix(STAR_ALIGN.out.versions.ifEmpty(null))

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
    // Split N Cigar Reads
    //

    SPLITNCIGAR (
        bam_markduplicates,
        fasta,
        fasta_fai,
        dict
    )

    bam_split = SPLITNCIGAR.out.bam

    ch_software_versions = ch_software_versions.mix(SPLITNCIGAR.out.versions.ifEmpty(null))


    //
    // Prepare Recalibration
    //

    PREPARE_RECALIBRATION (
        bam_split,
        params.use_gatk_spark,
        dict,
        fasta_fai,
        fasta,
        known_sites,
        known_sites_tbi,
        target_bed
    )

    table_bqsr = PREPARE_RECALIBRATION.out.table_bqsr
    bam_applybqsr = bam_split.join(table_bqsr)
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
        target_bed,
        gtf
    )

    bam_bqsr = RECALIBRATE.out.bam
    ch_software_versions = ch_software_versions.mix(RECALIBRATE.out.version.ifEmpty(null))

    //
    // HLATYPING
    //

    HLATYPING (
        bam_bqsr,
        hla_rna_fasta,
        input_sample
    )

    hla = HLATYPING.out.hla.groupTuple(by: 0).ifEmpty([])
    ch_software_versions = ch_software_versions.mix(HLATYPING.out.version.ifEmpty(null))

    //
    // VARIANT CALLING
    //

    RNA_VARIANT_CALLING (
        bam_bqsr,
        dbsnp,
        dbsnp_tbi,
        dict,
        fasta_fai,
        fasta
    )

    haplotypecaller_vcf = RNA_VARIANT_CALLING.out.haplotypecaller_vcf
    varscan_vcf         = RNA_VARIANT_CALLING.out.varscan_vcf

    ch_software_versions = ch_software_versions.mix(RNA_VARIANT_CALLING.out.version.ifEmpty(null))

    //
    // RNA COUNTS
    //

    FEATURECOUNTS (
        bam_bqsr,
        gtf
    )

    rna_counts = FEATURECOUNTS.out.counts
    rna_counts.map { meta, counts ->
        [meta.patient, counts]
    }.groupTuple(by: 0).set { rna_counts }

    ch_software_versions = ch_software_versions.mix(FEATURECOUNTS.out.version.ifEmpty(null))

    //
    // FILTERING
    //

    HAPLOTYPECALLER_FILTER (
        haplotypecaller_vcf,
        fasta,
        fasta_fai,
        dict
    )

    filtered_vcf = Channel.empty()

    filtered_vcf = filtered_vcf.mix(HAPLOTYPECALLER_FILTER.out.vcf)
    filtered_vcf = filtered_vcf.mix(varscan_vcf)

    filtered_vcf.groupTuple(by: 0).map { meta, paths ->
        [meta, *paths.sort( { it.getName().toString() })]
        }.set { vcf_to_merge }

    ch_software_versions = ch_software_versions.mix(HAPLOTYPECALLER_FILTER.out.version)

    //
    // COMBINEVARIANTS
    //

    COMBINE_VARIANTS (
        fasta,
        fasta_fai,
        dict,
        vcf_to_merge
    )

    merged_vcf = COMBINE_VARIANTS.out.vcf
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

    annotated_vcf = VEP.out.vcf.collect().flatten().collate(2).ifEmpty([])
    annotated_vcf.map{ meta, vcf ->
        meta.tumor = "${meta.patient}_${meta.sample}".toString()
        [meta.patient, meta.tumor, vcf]
    }.groupTuple(by: 0).set{ ann_vcf }

    ch_software_versions = ch_software_versions.mix(VEP.out.version.ifEmpty(null))

    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    GET_SOFTWARE_VERSIONS (
        input_sample,
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
        vep_vcf = ann_vcf
        counts  = rna_counts
        hla     = hla
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

        // TODO since it is mandatory: error/warning if not present?
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
