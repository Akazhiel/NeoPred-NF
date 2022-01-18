// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process YARA_INDEX {
    tag "$fasta"
    label 'process_medium'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'index', meta:[:], publish_by_meta:[]) }

    input:
    path fasta
    tuple val(meta), path(reads)


    output:
    path "yara"             , emit: index
    path "*.version.txt"    , emit: version

    script:
    def software = getSoftwareName(task.process)

    """
    mkdir yara
    yara_indexer \\
        $fasta \\
        -o "yara"
    mv *.{lf,rid,sa,txt}.* yara
    cp $fasta yara/yara.fasta

    echo \$(yara_indexer --version) | sed "s/^.*yara_indexer version: //; s/ .*\$//" > ${software}.version.txt
    """
}
