#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// params.outdir="qc_results"
// process qualityCheck{
//     // publishDir  "results", mode:'copy'
//     input:
//     tuple val(fastq), path(samples) 

//     output:
//     tuple val(fastq),path("results") 

//     script:
//     """
//     mkdir results
//     fastqc $fastq -o results
//     multiqc -f results -o results
//     """
// }

params.data="../Data/fastq_files"
params.outdir="qc_results"
path_ch=channel.fromFilePairs("${params.data}/*_{1,2}.fastq.gz")

params.out="../Data/ref"
path_ind_ch=channel.fromPath("../Data/ref/hs37d5.fa")

process qualityCheck{
    publishDir  "$params.outdir",pattern:"*.html", mode:'copy'
    input:
    path reads

    output:
    path 'results/*'

    script:
    """
    mkdir results
    fastqc -o results $reads
    multiqc -f results -o results
    """
}
fastq_pairs_ch = channel.fromPath("../Data/fastq_files/*.fastq.gz")



process filtering{
    input:
    tuple val(sample_id),path(fastq)
    output:
    tuple val(sample_id)
    script:
    """
    mkdir fastp_out
    fastp -i ${fastq[0]} -I ${fastq[1]} -o ${sample_id}_FP_R1.fq.gz -O ${sample_id}_FP_R2.fq.gz
    """
}
fastp_ch = channel.fromFilePairs("../Data/fastq_files/*.fastq.gz")
workflow{
    // qualityCheck(fastq_pairs_ch)
    filtering(path_ch)
}