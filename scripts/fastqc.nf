#!/usr/bin/env nextflow
params.data="../Data/fastq_files"
params.outdir="qc_results"
path_ch=channel.fromPath("${params.data}/*.fastq.gz")

params.out="../Data/ref"
path_ind_ch=channel.fromPath("../Data/ref/hs37d5.fa")

process qualityCheck{
    publishDir  "$params.outdir" mode:'copy'
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

// process qualityCheck{
//     publishDir  "$params.outdir" mode:'copy'
//     input:
//     file fastq from path_ch

//     output:
//     path("results") into qc_ch

//     script:
//     """
//     mkdir results
//     fastqc -o results $fastq
//     multiqc -f results -o results
//     """
// }


process filtering{
    input:
    output:
    script:
}

process index_ref{
    publishDir "$params.out"
    input:

    file ref from path_ind_ch
    output:
    path("../Data/ref/") into index_ch

    script:
    """
    bwa index $ref
    """
}