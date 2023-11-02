#!/usr/bin/env nextflow

// Processes to perform de novo assembly and annotate those assembled sequences
nextflow.enable.dsl=2
// Default parameters

params.min_coverage = 80 // linclust and reference genome alignment
// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.help = false
params.output = './results/'
params.alleles = null
params.allele_info = null

// Make sure that --output ends with trailing "/" characters
if (!params.output.endsWith("/")){
    params.output_folder = params.output.concat("/")
} else {
    params.output_folder = params.output
}

// Containers
container__anndata = "golob/python-anndata:0.9.2"


include { MMSeqs2_Cluster as MMSeqs2_Cluster_100 } from "./mmseqs" addParams(
    min_identity: 100,
    min_coverage: params.min_coverage
)
include { MMSeqs2_Cluster as MMSeqs2_Cluster_90 } from "./mmseqs" addParams(
    min_identity: 90,
    min_coverage: params.min_coverage
)
include { MMSeqs2_Cluster as MMSeqs2_Cluster_50 } from "./mmseqs" addParams(
    min_identity: 50,
    min_coverage: params.min_coverage
)

include { DiamondDB as DiamondIndex_C100 } from "./quantify" addParams(
    output_folder: params.output_folder
)

include { DiamondDB as DiamondIndex_C90 } from "./quantify" addParams(
    output_folder: params.output_folder
)

include { DiamondDB as DiamondIndex_C50 } from "./quantify" addParams(
    output_folder: params.output_folder
)

workflow Allele_clustering {
    take:
        alleles_faa_f
        allele_info_f

    main:

    
    // Cluster at 80 / 100 c/i
    MMSeqs2_Cluster_100(
        alleles_faa_f
    )

    MMSeqs2_Cluster_90(
        MMSeqs2_Cluster_100.out.seqs
    )

    MMSeqs2_Cluster_50(
        MMSeqs2_Cluster_90.out.seqs
    )
    // Summarize!

    SummarizeAllelesAndClusters(
        allele_info_f,
        MMSeqs2_Cluster_100.out.clusters,
        MMSeqs2_Cluster_90.out.clusters,
        MMSeqs2_Cluster_50.out.clusters,
        MMSeqs2_Cluster_100.out.seqs,
        MMSeqs2_Cluster_90.out.seqs,
        MMSeqs2_Cluster_50.out.seqs,        
    )

    DiamondIndex_C100(
        MMSeqs2_Cluster_100.out.seqs
    )

    DiamondIndex_C90(
        MMSeqs2_Cluster_90.out.seqs
    )

    DiamondIndex_C50(
        MMSeqs2_Cluster_50.out.seqs
    )

    emit:
        alleles = alleles_faa_f
        allele_info =  SummarizeAllelesAndClusters.out.AlleleClusterInfo
        centroids_C100 =  SummarizeAllelesAndClusters.out.C100
        dmdb_C100 = DiamondIndex_C100.out
        centroids_C90 = SummarizeAllelesAndClusters.out.C90
        dmdb_C90 = DiamondIndex_C90.out
        centroids_C50 = SummarizeAllelesAndClusters.out.C50
        dmdb_C50 = DiamondIndex_C50.out
    

    // */
}



process SummarizeAllelesAndClusters {
    tag "Summarize all of the avaliable alleles and output"
    container "${container__anndata}"
    label 'mem_medium'
    errorStrategy 'finish'
    publishDir "${params.output_folder}/alleles/", mode: "copy"
    
    input:
        path Alleles_csv
        path C100_tsv
        path C90_tsv
        path C50_tsv
        path Centroids_C100
        path Centroids_C90
        path Centroids_C50
    
    output:
        path('Allele_Cluster_info.csv.gz'), emit: AlleleClusterInfo
        path Centroids_C100, emit: C100
        path Centroids_C90, emit: C90
        path Centroids_C50, emit: C50        

"""
#!/usr/bin/env python3
import pandas as pd
a_i = pd.read_csv('${Alleles_csv}')
A_c100 = {
    r.allele: r.C100
    for i, r in 
    pd.read_csv('${C100_tsv}', delimiter='\\t', names=['C100', 'allele']).iterrows()
}
a_i['C100'] = a_i.allele.apply(A_c100.get)
c100_c90 = {
    r.C100: r.C90
    for i, r in 
    pd.read_csv('${C90_tsv}', delimiter='\\t', names=['C90', 'C100']).iterrows()
}
a_i['C90'] = a_i.C100.apply(c100_c90.get)
c90_c50 = {
    r.C90: r.C50
    for i, r in 
    pd.read_csv('${C50_tsv}', delimiter='\\t', names=['C50', 'C90']).iterrows()
}
a_i['C50'] = a_i.C90.apply(c90_c50.get)
a_i.to_csv(
    "Allele_Cluster_info.csv.gz",
    index=None
)
"""
}



//
// Steps to run gene-catalog independently.
//


// imports
include { Read_manifest } from './general'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run modules/allele_clustering.nf <ARGUMENTS>

    From a set of alleles and associated metadata for each allele (as below)
    generate clusters at 100%, 90% and 50% identity and associated diamond indicies

    Required Arguments:
      --alleles             faa.gz formatted allele sequences
      --allele_info         CSV formatted file with at least one column ('allele')
                                with at least one row per allele in the associated FAA

    Options:
      --output              Folder to place analysis outputs (default ./results)
      -w                    Working directory. Defaults to `./work`

    For Assembly:
      --min_coverage        For clustering. Default: 80 (same as uniref)
    """.stripIndent()
}


workflow {
    main:

 
    // Show help message if the user specifies the --help flag at runtime
    if (params.help || params.alleles == null || params.allele_info == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }

    Allele_clustering(
        file(params.alleles),
        file(params.allele_info)
  )

}
