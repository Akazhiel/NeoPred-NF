//
// MHC AFFINITY BINDING PREDICTION
//

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MHCFLURRY {
    tag "MHC_Binding_prediction"
    label 'process_low'
    publishDir "${params.outdir}/${patient}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'mhc_predict', meta:patient, publish_by_meta:[]) }

    input:
    tuple val(patient), path(hlas), path(variants)
    path(alleles)
    val(results)
    val(results_filter)
    val(seq_mode)
    val(cutoff)

    output:
    path("*_wt.csv"), optional:true, emit: wt_predictions
    path("*_mut.csv"), optional:true, emit: mut_predictions

    script:

    def software    = getSoftwareName(task.process)
    hla = "--hla " + hlas.join(' ')
    
    """
    mhc_predict.py \\
        $hla \\
        --variants $variants \\
        --alleles $alleles \\
        --mode $seq_mode \\
        --cutoff $cutoff  

    mhcflurry-predict-scan \\
        --alleles \$(cat allowed_alleles.txt) \\
        --results-${results} ${results_filter} \\
        --out predictions_wt.csv \\
         $options.args \\
        protein_sequences_wt.fasta  
    
    mhcflurry-predict-scan \\
        --alleles \$(cat allowed_alleles.txt) \\
        --results-${results} ${results_filter} \\
        --out predictions_mut.csv \\
         $options.args \\
        protein_sequences_mu.fasta 
    """
}