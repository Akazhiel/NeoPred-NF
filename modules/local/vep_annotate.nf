// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process VEP {
    tag "$meta.id"
    label 'process_high'
    // publishDir "${params.outdir}",
    //     mode: params.publish_dir_mode,
    //     saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['patient']) }

    container "hla-vep:latest"

    input:
    tuple val(meta), path(vcf), path(idx)
    path(fasta)
    path(cache)
    val(cache_version)
    val(vep_genome)

    output:
    tuple val(meta), path("*.ann.vcf")  , emit: vcf
    path "*.version.txt"                , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def dir_cache = cache ? "\${PWD}/${cache}" : "/.vep"

    """

    bcftools norm -m- $vcf | bcftools view -e 'ALT=="*"' -o ${meta.id}_split.vcf

    vep \\
        -i ${meta.id}_split.vcf \\
        -o ${prefix}.ann.vcf \\
        $options.args \\
        --assembly ${params.genome} \\
        --cache \\
        --verbose \\
        --fasta ${fasta} \\
        --cache_version $cache_version \\
        --dir_cache $dir_cache \\
        --fork $task.cpus \\
        $options.args2

    echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//' > ${software}.version.txt
    """
}
