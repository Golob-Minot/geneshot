#!/usr/bin/env nextflow

nextflow.enable.dsl=2

container__mmseqs2 = 'golob/mmseqs:20230807A'

process MMSeqs2_Cluster {
    tag "Use mmseqs2 to cluster. Non-linear"
    container = container__mmseqs2 
    label 'mem_veryhigh'
    errorStrategy 'ignore'
    // publishDir "${params.output_folder}/cluster/", mode: "copy"



    input:
        path "original.fastp.gz"

    output:
        path "clusters.C${params.min_identity}.tsv.gz", emit: clusters
        path "centroids.C${params.min_identity}.fastp.gz", emit: seqs

"""
set -e

mkdir db/
mkdir cluster/
entrypoint createdb original.fastp.gz db/original.DB
entrypoint  cluster db/original.DB cluster/cluster.DB /tmp/mmseq -c ${params.min_coverage / 100}  --min-seq-id ${params.min_identity / 100}
entrypoint  createtsv db/original.DB db/original.DB cluster/cluster.DB clusters.C${params.min_identity}.tsv
entrypoint  createsubdb cluster/cluster.DB db/original.DB cluster/cluster_rep.DB
entrypoint  convert2fasta cluster/cluster_rep.DB centroids.C${params.min_identity}.fastp
gzip centroids.C${params.min_identity}.fastp
gzip clusters.C${params.min_identity}.tsv
"""

}



//
// Steps to run mmseqs2 independently.
//

params.help = false
params.in_fastp_gz = null
params.output = '.'
params.min_coverage = 80
params.min_identity = 50

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run Golob-Minot/geneshot/modules/mmseqs.nf <ARGUMENTS>
    
    Required Arguments:
      --in_fastp_gz         gzipped FASTP format file to cluster
      --output              Directory into which outputs should be placed (default .)

    Options:
      --min_coverage        Length cutoff used to combine similar genes (default: 80)
      --min_identity        Percent identity (default: 50)
    """.stripIndent()
}

workflow {
    main:

 
    // Show help message if the user specifies the --help flag at runtime
    if (params.help || params.in_fastp_gz == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }
    MMSeqs2_Cluster(
        file(params.in_fastp_gz)
    )



}
