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

dna_tumor_cov      = Channel.value(params.dna_tumor_cov)
dna_tumor_depth    = Channel.value(params.dna_tumor_depth)
dna_tumor_vaf      = Channel.value(params.dna_tumor_vaf)
dna_normal_cov     = Channel.value(params.dna_normal_cov)
dna_normal_vaf     = Channel.value(params.dna_normal_vaf)
tumor_normal_ratio = Channel.value(params.tumor_normal_ratio)
dna_snv_callers    = Channel.value(params.dna_snv_callers)
dna_indel_callers  = Channel.value(params.dna_indel_callers)
rna_tumor_cov      = Channel.value(params.rna_tumor_cov)
rna_tumor_depth    = Channel.value(params.rna_tumor_depth)
rna_tumor_vaf      = Channel.value(params.rna_tumor_vaf)
rna_callers        = Channel.value(params.rna_callers)
cdna_dict          = Channel.fromPath(params.cdna_dict).collect()
aa_dict            = Channel.fromPath(params.aa_dict).collect()
pyensembl_version  = Channel.value(params.pyensembl)

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

include { MERGE_VARIANTS }        from '../modules/local/merge_variants'        addParams(options: [:])

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
    vep_vcf_dna
    vep_vcf_rna
    rna_counts

    main:

    if (vep_vcf_dna && vep_vcf_rna) {
        vcfs = vep_vcf_dna.join(vep_vcf_rna, remainder: true)
        vcfs = vcfs.join(rna_counts, remainder: true).groupTuple(by: 0, size: 6, remainder: true)

        vcfs.branch{
            vcfs_5: it.size() == 5
            vcfs_6: it.size() == 6
        }.set { vcfs }

        vcfs_5 = vcfs.vcfs_5.map { patient, dna, dna_vcf, rna, rna_vcf ->
            [patient, *dna, *dna_vcf, [], [], []]
            }
        
        vcfs_6 = vcfs.vcfs_6.map { patient, dna, dna_vcf, rna, rna_vcf, rna_counts ->
            [patient, *dna, *dna_vcf, *rna, *rna_vcf, *rna_counts]
            }

        vcfs_to_merge = vcfs_5.mix(vcfs_6)
    } else if (vep_vcf_dna && !vep_vcf_rna) {
        vcfs_to_merge = vep_vcf_dna.map { patient, dna, dna_vcf -> 
        [patient, dna, dna_vcf, [], [], []]
        }
    } else if (!vep_vcf_dna && vep_vcf_rna) {
        vcfs_to_merge = vep_vcf_rna.join(rna_counts).map { patient, rna, rna_vcf, rna_counts -> 
        [patient, [], [], rna, rna_vcf, rna_counts]
        }
    }

    MERGE_VARIANTS (
        vcfs_to_merge,
        dna_tumor_cov,     
        dna_tumor_depth,   
        dna_tumor_vaf,     
        dna_normal_cov,    
        dna_normal_vaf,    
        tumor_normal_ratio,
        dna_snv_callers,   
        dna_indel_callers, 
        rna_tumor_cov,     
        rna_tumor_depth,   
        rna_tumor_vaf,     
        rna_callers,
        cdna_dict,
        aa_dict    
    )

    emit:

    overlap_final         = MERGE_VARIANTS.out.overlap_final
    overlap_discarded     = MERGE_VARIANTS.out.overlap_discarded
    overlap_rna           = MERGE_VARIANTS.out.overlap_rna
    overlap_rna_discarded = MERGE_VARIANTS.out.overlap_rna_discarded  
    final_counts          = MERGE_VARIANTS.out.final_counts
}