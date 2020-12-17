#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

  This subworkflow runs Beta-binomial Distribution Modeling (BBDML) on the CAGs.
  You need to provide:
  - results.hdf5
  - directory containing the abund/ subfolder
  - A matrix of exogenous covariates for both abundance AND dispersion
*/

// Using DSL-2
nextflow.preview.dsl=2

// Parameters
params.input_hdf = false
params.input_folder = false
params.X = false
params.X_star = false
params.output_folder = false
params.output_prefix = ""
params.batchsize = false
params.help = false


params.fdr_method = "fdr_bh"



// Docker containers
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run corncob with a specified formula on a geneshot output file.
    This utility will create a new geneshot output file in which
    the only difference is the corncob output results stored within.

    Usage:

    nextflow run Golob-Minot/geneshot/run_bbdml.nf <ARGUMENTS>
    
    Options:
      --input_hdf           Geneshot results HDF5 file to be used for analysis
      --input_folder        Folder containing geneshot output to be used for analysis (must contain the subfolder "abund/")
      --output_folder       Folder where the regression results should go
      --output_prefix       Prefix to prepend to outputs. Default is empty ""
      --batchsize           OPTIONAL. Break CAGs into batches of this size for parallel analysis
      -w                    Working directory. Defaults to `./work`
      
    For Statistical Analysis:
      --X                   Exogenous covariates for abundance
      --X_star              Exogenous covariates for dispersion
      --fdr_method          FDR method used to calculate q-values for associations (default: 'fdr_bh')
      

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.input_hdf == false || params.input_folder == false || params.output_folder == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Import the process used to aggregate readcounts across all samples
// This process is only run if `abund/CAG.readcounts.T.csv.gz` can not be found yet
// After the process completes, a new `abund/CAG.readcounts.T.csv.gz` will be added to that input folder
include { ExtractCountsT } from './modules/statistics' params(
    output_folder: params.input_folder
)

include { 
    RunBBDML
    RunBBDML_nodisp
    RunBBDML_noabund
    Shard_CAG_Readcounts
    Join_BBDML
    Add_FDR
} from './modules/statistics' params(
    output_folder: params.output_folder,
    output_prefix: params.output_prefix,
    batchsize: params.batchsize
)
// *** //



workflow {

    // Make sure we can find the input file
    if(file(params.input_hdf).isEmpty()){
        log.info"""Cannot find input file ${params.input_hdf}""".stripIndent()
        exit 0
    }

    // Check to see if `abund/CAG.readcounts.csv.gz` is present in the input folder
    // For context, this file is only created by geneshot for runs which originally
    // included a `--formula`. 
    // If this file cannot be found, then we will run a process to create that file
    // and also to publish that file back to the input folder.
    abund_folder = "${params.input_folder.replaceAll(/\/$/, "")}/abund"
    if(file("${abund_folder}/CAG.readcounts.T.csv.gz").isEmpty()){
        extractCountsT(
            Channel.fromPath(
                "${abund_folder}/details/*json.gz"
            ).toSortedList(),
            file(params.input_hdf)
        )
        crc_csv = extractCountsT.out
    } else {
        crc_csv = file("${abund_folder}/CAG.readcounts.T.csv.gz")
    }

    // Shard or chunk CAGs
    if (params.batchsize == false) {
        crc_blocks = Channel.from(
            crc_csv
        )

    } else {
        Shard_CAG_Readcounts(
            crc_csv
        )
        crc_blocks = Shard_CAG_Readcounts.out.flatten()
    }

    // Run BBDML on each chunk
    if (params.X == false && params.X_star == false){
        log.info("Hmmm. No covariates (X or X_star) were provided. Are you sure this is what you wanted?")
        return None
    } else if (params.X_star == false) {
        RunBBDML_nodisp(
            crc_blocks,
            file(params.X),
        )
        bbdml_out = RunBBDML_nodisp.out
    } else if (params.X == false) {
        RunBBDML_noabund(
            crc_blocks,
            file(params.X_star)
        )    
        bbdml_out = RunBBDML_noabund.out
    } else {
        RunBBDML(
            crc_blocks,
            file(params.X),
            file(params.X_star)
        )
        bbdml_out = RunBBDML.out

    }

    Join_BBDML(
        bbdml_out.toSortedList()
    )

    Add_FDR(
        Join_BBDML.out
    )

    emit:
        Add_FDR.out
}