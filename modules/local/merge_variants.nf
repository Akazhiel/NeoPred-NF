//
// MERGE DNA AND RNA VARIANTS
//

include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process MERGE_VARIANTS {
    tag "Merge_variants"
    label 'process_low'
    publishDir "${params.outdir}/${patient}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'overlap_merge', meta:patient, publish_by_meta:[]) }

    input:
    tuple val(patient), val(tumor_dna), path(vcfs_dna), val(tumor_rna), path(vcfs_rna), path(counts)
    val dna_tumor_cov     
    val dna_tumor_depth   
    val dna_tumor_vaf     
    val dna_normal_cov    
    val dna_normal_vaf    
    val tumor_normal_ratio
    val dna_snv_callers   
    val dna_indel_callers 
    val rna_tumor_cov     
    val rna_tumor_depth   
    val rna_tumor_vaf     
    val rna_callers    
    path cdna_dict
    path aa_dict   

    output:
    tuple val(patient), path("*_final.txt"),                        emit: overlap_final
    tuple val(patient), path("*_final_rna_unique.txt"),             emit: overlap_rna
    tuple val(patient), path("*_final_discarded.txt"),              emit: overlap_discarded
    tuple val(patient), path("*_rna_unique_discarded.txt"),         emit: overlap_rna_discarded
    tuple val(patient), path("*.gene.counts.final"), optional:true, emit: final_counts

    script:
    def software    = getSoftwareName(task.process)
    def vcf_dna     = vcfs_dna  ? "--dna " + vcfs_dna.join(' ')         : ''
    def samples_dna = tumor_dna ? "--dna-names " + tumor_dna.join(' ')  : ''
    def vcf_rna     = vcfs_rna  ? "--rna " + vcfs_rna.join(' ')         : ''
    def samples_rna = tumor_rna ? "--rna-names " + tumor_rna.join(' ')  : ''
    def rna_counts  = counts    ? "--rna-counts " + counts.join(' ')    : ''

    """
    merge_variants.py \\
        $vcf_dna \\
        $samples_dna \\
        $vcf_rna \\
        $samples_rna \\
        $rna_counts \\
        --filter-dna-tumor-cov $dna_tumor_cov \\
        --filter-dna-tumor-depth $dna_tumor_depth \\
        --filter-dna-tumor-vaf $dna_tumor_vaf \\
        --filter-dna-normal-cov $dna_normal_cov \\
        --filter-dna-normal-vaf $dna_normal_vaf \\
        --filter-dna-tn-ratio $tumor_normal_ratio \\
        --filter-dna-snv-callers $dna_snv_callers \\
        --filter-dna-indel-callers $dna_indel_callers \\
        --filter-rna-tumor-cov $rna_tumor_cov \\
        --filter-rna-tumor-depth $rna_tumor_depth \\
        --filter-rna-tumor-vaf $rna_tumor_vaf \\
        --filter-rna-callers $rna_callers \\
        --ensembl-version $params.pyensembl \\
        --dictAA $aa_dict \\
        --dictcDNA $cdna_dict
    """
}