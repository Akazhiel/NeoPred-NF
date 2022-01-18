#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)
params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')
params.fasta_fai = WorkflowMain.getGenomeAttribute(params, 'fasta_fai')
params.star_index = WorkflowMain.getGenomeAttribute(params, 'star_index')
params.gtf = WorkflowMain.getGenomeAttribute(params, 'gtf')
params.dict = WorkflowMain.getGenomeAttribute(params, 'dict')
params.bwa = WorkflowMain.getGenomeAttribute(params, 'bwa')
params.germline = WorkflowMain.getGenomeAttribute(params, 'germline')
params.germline_index = WorkflowMain.getGenomeAttribute(params, 'germline_index')
params.pon = WorkflowMain.getGenomeAttribute(params, 'pon')
params.pon_index = WorkflowMain.getGenomeAttribute(params, 'pon_index')
params.dbsnp = WorkflowMain.getGenomeAttribute(params, 'dbsnp')
params.dbsnp_index = WorkflowMain.getGenomeAttribute(params, 'dbsnp_index')
params.intervals = WorkflowMain.getGenomeAttribute(params, 'intervals')
params.vep_genome = WorkflowMain.getGenomeAttribute(params, 'vep_genome')
params.vep_cache_version = WorkflowMain.getGenomeAttribute(params, 'vep_cache_version')
params.known_indels = WorkflowMain.getGenomeAttribute(params, 'known_indels')
params.known_indels_index = WorkflowMain.getGenomeAttribute(params, 'known_indels_index')
params.hla_reference_dna = WorkflowMain.getGenomeAttribute(params, 'hla_reference_dna')
params.hla_reference_rna = WorkflowMain.getGenomeAttribute(params, 'hla_reference_rna')
params.aa_dict = WorkflowMain.getGenomeAttribute(params, 'AA_dict')
params.cdna_dict = WorkflowMain.getGenomeAttribute(params, 'cDNA_dict')
params.pyensembl = WorkflowMain.getGenomeAttribute(params, 'pyensembl')
params.alleles = WorkflowMain.getGenomeAttribute(params, 'alleles')

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMain.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { NEOPRED_DNA }     from './workflows/neopreddna'
include { NEOPRED_RNA }     from './workflows/neopredrna'
include { MERGE_RESULTS }   from './workflows/merge_results'
include { MHC_PREDICT }     from './workflows/mhc_predict'

//
// WORKFLOW: Run main nf-core/neoprednf analysis pipeline
//
workflow NFCORE_NEOPREDNF {

    if ( params.DNA && params.RNA ) {
        NEOPRED_DNA ()
        NEOPRED_RNA ()
        MERGE_RESULTS (
            NEOPRED_DNA.out.vep_vcf,
            NEOPRED_RNA.out.vep_vcf,
            NEOPRED_RNA.out.counts
        )
        MHC_PREDICT (
            NEOPRED_DNA.out.hla,
            NEOPRED_RNA.out.hla,
            MERGE_RESULTS.out.overlap_final
        )
    } else if ( params.DNA && !params.RNA) {
        NEOPRED_DNA ()
        MERGE_RESULTS (
            NEOPRED_DNA.out.vep_vcf,
            [],
            []
        )
        MHC_PREDICT (
            NEOPRED_DNA.out.hla,
            [],
            MERGE_RESULTS.out.overlap_final
        )
    } else if ( !params.DNA && params.RNA) {
        NEOPRED_RNA ()
        MERGE_RESULTS (
            [],
            NEOPRED_RNA.out.vep_vcf,
            NEOPRED_RNA.out.counts
        )
        MHC_PREDICT (
            [],
            NEOPRED_RNA.out.hla,
            MERGE_RESULTS.out.overlap_rna
        )
    }

}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_NEOPREDNF ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
