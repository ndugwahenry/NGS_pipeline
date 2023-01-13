#!/usr/bin/env nextflow
params.data="../Data/fastq_files"
params.outdir="qc_results"
path_ch=channel.fromPath("${params.data}/*.fastq.gz")

params.out="../Data/ref"
path_ind_ch=channel.fromPath("../Data/ref/hs37d5.fa")

process qualityCheck{
    publishDir "$params.outdir"
    input:
    file fastq from path_ch

    output:
    path("results") into qc_ch

    script:
    """
    mkdir results
    fastqc -o results $fastq
    multiqc -f results -o results
    """
}

process index_ref{
    publishDir "$params.out"
    input:

    file ref from path_ind_ch
    output:
    path("../Data/ref/")

    script:
    """
    bwa index $ref
    """
}