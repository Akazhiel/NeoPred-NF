#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/neoprednf
========================================================================================
    Github : https://github.com/nf-core/neoprednf
    Website: https://nf-co.re/neoprednf
    Slack  : https://nfcore.slack.com/channels/neoprednf
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/

params.fasta = WorkflowMain.getGenomeAttribute(params, 'fasta')

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

include { NEOPREDNF } from './workflows/neoprednf'

//
// WORKFLOW: Run main nf-core/neoprednf analysis pipeline
//
workflow NFCORE_NEOPREDNF {
    NEOPREDNF ()
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
