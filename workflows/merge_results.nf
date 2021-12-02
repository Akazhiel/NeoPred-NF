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

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()

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

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
// def multiqc_report = []

workflow MERGE_RESULTS {
    take:
    vep_vcf

    // vep_vcf.view()

}
