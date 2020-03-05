#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

    This is a workflow oriented around obtaining the composition of a community via WGS data.
    Reuses components (primarily pre-processing) from the broader geneshot. 
*/

// Using DSL-2
nextflow.preview.dsl=2

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.nopreprocess = false
params.savereads = false
params.help = false
params.output = './results'
params.output_prefix = 'geneshot'
params.manifest = null

// Preprocessing options
params.adapter_F = "CTGTCTCTTATACACATCT"
params.adapter_R = "CTGTCTCTTATACACATCT"
params.hg_index_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30


// Function which prints help message text
def helpMessage() {
    log.info"""
    A workflow oriented around obtaining the composition of a community via WGS data.
    Reuses components (primarily pre-processing) from the broader geneshot. 

    Usage:

    nextflow run composition_only.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    Options:
      --output              Folder to place analysis outputs (default ./results)
      --output_prefix       Text used as a prefix for summary HDF5 output files (default: geneshot)
      --nopreprocess        If specified, omit the preprocessing steps (removing adapters and human sequences)
      --savereads           If specified, save the preprocessed reads to the output folder (inside qc/)
      -w                    Working directory. Defaults to `./work`

    For preprocessing:
      --hg_index_url        URL for human genome index, defaults to current HG
      --hg_index            Cached copy of the bwa indexed human genome, TGZ format
      --adapter_F           Forward sequencing adapter sequence (to be removed)
      --adapter_R           Reverse sequencing adapter sequence (to be removed)
                              (Adapter sequences default to nextera adapters)
      --min_hg_align_score  Minimum alignment score for human genome (default 30)
    
    Manifest file:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This can be repeated. 
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If index reads are provided, the column titles should be 'I1' and 'I2'

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.manifest == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Make sure that --output ends with trailing "/" characters
if (!params.output.endsWith("/")){
    output_folder = params.output.concat("/")
} else {
    output_folder = params.output
}

// Import the preprocess_wf module
include read_manifest from './modules/preprocess'
include preprocess_wf from './modules/preprocess' params(
    manifest: params.manifest,
    adapter_F: params.adapter_F,
    adapter_R: params.adapter_R,
    hg_index: params.hg_index,
    hg_index_url: params.hg_index_url,
    min_hg_align_score: params.min_hg_align_score,
)
// Import some general tasks, such as combineReads and writeManifest
include countReads from './modules/general'
include countReadsSummary from './modules/general' params(
    output_folder: output_folder
)

include writeManifest from './modules/general' params(
    savereads: params.savereads,
    output_folder: output_folder
)
include combineReads from './modules/general' params(
    savereads: params.savereads,
    output_folder: output_folder
)


workflow {
    main:

    // Phase 0: Validation of input data

    // If the user specifies a `--formula`, the first step in the process
    // will be to ensure that the formula is written correctly, and is
    // compatible with the data provided in the manifest
    if ( params.formula ) {
        validation_wf(
            file(params.manifest)
        )
        manifest_file = validation_wf.out
    } else {
        manifest_file = Channel.from(file(params.manifest))
    }

    // Phase I: Preprocessing
    if (!params.nopreprocess) {

        // Run the entire preprocessing workflow
        preprocess_wf(
            manifest_file
        )

        // Combine the reads by specimen name
        combineReads(preprocess_wf.out.groupTuple())

    } else {
        // If the user specified --nopreprocess, then just 
        // read the manifest and combine by specimen
        combineReads(
            read_manifest(
                manifest_file
            ).filter { r ->
                (r.specimen != null) &&
                (r.R1 != null) &&
                (r.R2 != null) &&
                (r.specimen != "") &&
                (r.R1 != "") &&
                (r.R2 != "")
            }.filter {
                r -> (!file(r.R1).isEmpty() && !file(r.R2).isEmpty())
            }.map { 
                r -> [r.specimen, file(r.R1), file(r.R2)]
            }.groupTuple()
        )

    }

    // If the user specified --savereads, write out the manifest
    if (params.savereads) {
        writeManifest(
            combineReads.out
        )
    }

    // Count the reads for every sample individually (just take the first of the pair of reads)
    countReads(
        combineReads.out.map {
            r -> [r[0], r[1], r[2]]
        }
    )

    // Make a summary of every sample and write it out to --output
    countReadsSummary(
        countReads.out.collect()
    )

    // ##############################
    // # Composition via MetaPhlAn2 #
    // ##############################



    // ###################
    // # GATHER RESULTS #
    // ###################

    // Start by gathering all of the results which are generated
    // no matter what options were selected by the user
    // NOTE: The code used here is imported from ./modules/general.nf

    collectAbundances(
        countReadsSummary.out,
        manifest_file
    )

    // "Repack" the HDF5, which enhances space efficiency and adds GZIP compression
    repackFullHDF(
        finalHDF
    )

    publish:
        corncob_results to: "${output_folder}/stats/", enabled: params.formula
        alignment_wf.out.famli_json_list to: "${output_folder}/abund/details/"
        repackFullHDF.out to: "${output_folder}", mode: "copy", overwrite: true
        repackSummaryHDF.out to: "${output_folder}", mode: "copy", overwrite: true

}
