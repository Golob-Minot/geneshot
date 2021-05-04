#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

  This utility performs taxonomic and/or functional annotations on an existing geneshot
  results object
*/

// Using DSL-2
nextflow.enable.dsl=2

// Parameters
params.input_hdf = false
params.gene_fasta = false
params.output_hdf = false
params.output_folder = false
params.help = false
params.min_coverage = 50
params.min_identity = 90
params.taxonomic_dmnd = false
params.ncbi_taxdump = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
params.eggnog_db = false
params.eggnog_dmnd = false
params.noannot = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Add functional and/or taxonomic annotations to a geneshot output file.
    This utility will create a new geneshot output file in which
    the only difference is the annotation output results stored within.

    Usage:

    nextflow run Golob-Minot/geneshot/run_annotations.nf <ARGUMENTS>
    
    Options:
      --input_hdf           Location for input HDF
      --gene_fasta          Location for input 'genes.fasta.gz'
      --output_folder       Location for output HDF
      --output_hdf          Name of output HDF to be placed in the output folder
      --taxonomic_dmnd      Database used for taxonomic annotation (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd)
      --ncbi_taxdump        Reference describing the NCBI Taxonomy
                            (default: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
      --eggnog_dmnd         One of two databases used for functional annotation with eggNOG (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-06-17-eggNOG-v5.0/eggnog_proteins.dmnd)
      --eggnog_db           One of two databases used for functional annotation with eggNOG (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-06-17-eggNOG-v5.0/eggnog.db)
    
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.input_hdf == false || params.gene_fasta == false || params.output_hdf == false || params.output_folder == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Show help message if the user does not specify any annotations
if (params.taxonomic_dmnd == false && params.eggnog_dmnd == false && params.eggnog_db == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Import the process used to add annotation results to the output
include { readTaxonomy } from './modules/general'
include { addTaxResults } from './modules/general'

include { addEggnogResults } from './modules/general'

include { repackHDF } from './modules/general' params(
    output_folder: params.output_folder
)

// Import the workflows used for annotation
include { annotation_wf } from './modules/genecatalog' params(
    output_folder: params.output_folder,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    eggnog_db: params.eggnog_db,
    eggnog_dmnd: params.eggnog_dmnd,
    taxonomic_dmnd: params.taxonomic_dmnd,
    noannot: params.noannot
)

// Process to rename the input HDF file
process renameHDF{
    container "ubuntu:20.04"
    label 'io_limited'

    input:
        path "INPUT.RESULTS.HDF"
        val new_file_name

    output:
        path "${new_file_name}"

"""#!/bin/bash

mv INPUT.RESULTS.HDF ${new_file_name}
"""
}


workflow {

    // Make sure we can find the input files
    if(file(params.input_hdf).isEmpty()){
        log.info"""Cannot find input file ${params.input_hdf}""".stripIndent()
        exit 0
    }
    if(file(params.gene_fasta).isEmpty()){
        log.info"""Cannot find input file ${params.gene_fasta}""".stripIndent()
        exit 0
    }

    // Run the annotation steps on the gene catalog
    annotation_wf(
        file(params.gene_fasta)
    )

    // Rename the input HDF
    renameHDF(
        file(params.input_hdf),
        params.output_hdf
    )
    outputHDF = renameHDF.out

    // Add the eggNOG results
    if (annotation_wf.out.eggnog_tsv != false && params.eggnog_db != false && params.eggnog_dmnd != false) {
        addEggnogResults(
            outputHDF,
            annotation_wf.out.eggnog_tsv
        )
        outputHDF = addEggnogResults.out
    }
    
    // Add the taxonomic assignment results
    if (annotation_wf.out.tax_tsv != false && params.ncbi_taxdump != false && params.taxonomic_dmnd != false) {
        readTaxonomy(
            file(params.ncbi_taxdump)
        )

        addTaxResults(
            outputHDF,
            annotation_wf.out.tax_tsv,
            readTaxonomy.out
        )

        outputHDF = addTaxResults.out
    }

    // Repack the HDF
    repackHDF(
        outputHDF
    )

}