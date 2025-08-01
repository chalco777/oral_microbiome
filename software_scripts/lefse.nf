#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.fastqDir = "$launchDir/data/"

process Lefse {
    container "quay.io/biocontainers/lefse:1.0.7--1"
    publishDir 'results/lefse', mode: 'copy'

    input:
    tuple val(sample_id), path(query_file)

    output:
    val sample_id

    script:
    """
    lefse_format_input.py ${query_file} species_reads_tabla.in -c 1 -s -1 -u 2 -o 1000000

    lefse_run.py species_reads_tabla.in species_reads_tabla_L_2.0.res

    lefse_plot_res.py species_reads_tabla_L_2.0.res algae_lefse_LDA_2.0.png --dpi 1000

    """
}

workflow {
    channel_fastq = channel.fromFilePairs("$params.fastqDir/kracken_report*.tsv",size:1)
    //channel_fastq.view()
    Lefse(channel_fastq)
}
