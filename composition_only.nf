#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

    This is a workflow oriented around obtaining the composition of a community via WGS data.
    Reuses components (primarily pre-processing) from the broader geneshot. 
*/

// Using DSL-2
nextflow.enable.dsl=2

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.nopreprocess = false
params.savereads = false
params.help = false
params.output = './results'
params.manifest = false

// Preprocessing options
params.hg_index_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30

// Make sure that --output ends with trailing "/" characters
if (!params.output.endsWith("/")){
    output_folder = params.output.concat("/")
} else {
    output_folder = params.output
}


// Import the preprocess_wf module
include { Read_manifest } from './modules/general'
include { Preprocess_wf } from './modules/preprocess' params(
    hg_index: params.hg_index,
    hg_index_url: params.hg_index_url,
    min_hg_align_score: params.min_hg_align_score,
    output: output_folder,
    savereads: params.savereads
)
// Import some general tasks
include { CombineReads } from './modules/preprocess' params(
    savereads: params.savereads,
    output_folder: output_folder
)


// Import from composition_wf module
include { composition_wf } from './modules/composition' params(
    output_folder: output_folder
)

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
      --nopreprocess        If specified, omit the preprocessing steps (removing adapters and human sequences)
      --savereads           If specified, save the preprocessed reads to the output folder (inside qc/)
      -w                    Working directory. Defaults to `./work`

    For preprocessing:
      --hg_index_url        URL for human genome index, defaults to current HG
      --hg_index            Cached copy of the bwa indexed human genome, TGZ format
      --min_hg_align_score  Minimum alignment score for human genome (default 30)
    
    Manifest file:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This can be repeated. 
      Data for preprocessing is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If index reads are provided, the column titles should be `I1` and `I2`
      If you wish to provide already processed data in fasta format, please include it in `R1` alone,
        with only *one* file specified per specimen.
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || !params.manifest){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}


workflow {
    main:

    // Phase 0: Validation of input data
    manifest_file = Channel.from(file(params.manifest))
    // Read manifest splits out our manifest.
    manifest_qced = Read_manifest(manifest_file)

    // Phase I: Preprocessing
    if (!params.nopreprocess) {
        // Run the entire preprocessing workflow
        Preprocess_wf(
            manifest_qced.valid_paired_indexed,
             manifest_qced.valid_paired
        )
        // Combine the reads by specimen name
        CombineReads(Preprocess_wf.out.groupTuple())

    } else {
        // If the user specified --nopreprocess, then just 
        // read the manifest and combine by specimen
        CombineReads(
            manifest_qced.valid_paired.mix(manifest_qced.valid_paired_indexed)
            .map { 
                r -> [r.specimen, file(r.R1), file(r.R2)]
            }.groupTuple()
        )

    }

    // ################
    // # Composition  #
    // ################
    composition_wf(
        CombineReads.out,
        manifest_qced.valid_unpaired.map{ r-> 
            [r.specimen, file(r.R1)]
        }
    )

}
