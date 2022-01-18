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

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

ch_dummy_file = Channel.fromPath("$projectDir/assets/dummy_file.txt", checkIfExists: true).collect()

alleles         = Channel.fromPath(params.alleles).collect()
results         = Channel.value(params.results)
results_filter  = Channel.value(params.results_filter)
cutoff          = Channel.value(params.cutoff)
seq_mode        = Channel.value(params.seq_mode)

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

include { MHCFLURRY }             from '../modules/local/mhcflurry'             addParams( options: modules['mhcflurry'] )

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

workflow MHC_PREDICT {
    take:
    hla_dna
    hla_rna
    variants

    main:

    if (hla_dna && hla_rna) {
        hla = hla_dna.join(hla_rna, remainder: true).join(variants, remainder: true)

        hla.branch {
            not_null: it[2] != null
            is_null: it[2] == null
        }.set { hlas }
        
        hla_not_null = hlas.not_null.map { patient, hla_dna, hla_rna, variants ->
            [patient, [*hla_dna, *hla_rna], variants]    
        }.dump()

        hla_null = hlas.is_null.map { patient, hla_dna, hla_rna, variants ->
            [patient, hla_dna, variants]    
        }

        hla_mhc = hla_not_null.mix(hla_null)
               
    } else if (hla_dna && !hla_rna) {
        hla_mhc = hla_dna.join(variants).map { patient, hla_dna, variants -> 
        [patient, hla_dna, variants]
        }
    } else if (!hla_dna && hla_rna) {
        hla_mhc = hla_rna.join(variants).map { patient, hla_rna, variants -> 
        [patient, hla_rna, variants]
        }
    }

    MHCFLURRY (
        hla_mhc,
        alleles,
        results,
        results_filter,
        seq_mode,
        cutoff
    )
}