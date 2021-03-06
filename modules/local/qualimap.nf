// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process QUALIMAP_BAMQC {
    tag "$meta.id"
    label 'process_medium'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:[:]) }

    input:
    tuple val(meta), path(bam), path(bai)
    path gff
    val use_gff

    output:
    tuple val(meta), path("${prefix}"), emit: results
    path  "*.version.txt"             , emit: version

    script:
    def software       = getSoftwareName(task.process)
    prefix             = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    def collect_pairs = meta.single_end ? '' : '--collect-overlap-pairs'
    def memory        = task.memory.toGiga() + "G"
    def regions       = use_gff ? "--gff $gff" : ''

    def strandedness = 'non-strand-specific'

    """
    unset DISPLAY
    mkdir tmp
    export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
    qualimap \\
        --java-mem-size=$memory \\
        bamqc \\
        $options.args \\
        -bam $bam \\
        $regions \\
        -p $strandedness \\
        $collect_pairs \\
        -outdir $prefix \\
        -nt $task.cpus
    echo \$(qualimap 2>&1) | sed 's/^.*QualiMap v.//; s/Built.*\$//' > ${software}.version.txt
    """
}
