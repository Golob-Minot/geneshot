#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

  There is an optional preprocessing submodule, 
  followed by assembly (spades) + extraction of peptide coding sequences from the contigs
  The short reads are then aligned against the assembled peptides plus uniref100.
  We use the FAMLI algorithm to adjuticate these alignments.
  Annotations can follow. 

  I. (Optional) Geneshot preprocessing submodule:
    Steps:
    1) (if index is available): barcodecop to verify demultiplexing
    2) cutadapt to remove adapters.
    3) remove human reads via
      3A) downloading the cached human genome index
      3B) aligning against the human genome and extracting unpaired reads
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
params.manifest = null

// Preprocessing options
params.adapter_F = "CTGTCTCTTATACACATCT"
params.adapter_R = "CTGTCTCTTATACACATCT"
params.hg_index_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30

// Assembly options
params.phred_offset = 33 // spades
params.centre = 'geneshot' // prokka
params.min_identity = 90 // mmseqs2
params.min_coverage = 50 // mmseqs2
params.dmnd_min_identity = 80 // DIAMOND
params.dmnd_min_coverage = 50 // DIAMOND
params.dmnd_top_pct = 1 // DIAMOND
params.dmnd_min_score = 20 // DIAMOND
params.gencode = 11 //DIAMOND

// Annotation options
params.noannot = false
params.taxonomic_dmnd = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd"
params.eggnog_dmnd = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog.db"
params.eggnog_db = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog_proteins.dmnd"


// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/geneshot <ARGUMENTS>
    
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
      --adapter_F           Forward sequencing adapter sequence (to be removed)
      --adapter_R           Reverse sequencing adapter sequence (to be removed)
                              (Adapter sequences default to nextera adapters)
      --min_hg_align_score  Minimum alignment score for human genome (default 30)
    For Assembly:
      --phred_offset        for spades. Default 33.
      --centre              Centre for use in prokka. default = 'geneshot'
      --min_identity        Amino acid identity cutoff used to combine similar genes (default: 90) (mmseqs2)
      --min_coverage        Length cutoff used to combine similar genes (default: 50) (mmseqs2)
      --dmnd_min_identity   Amino acid identity cutoff used to align short reads (default: 90) (DIAMOND)
      --dmnd_min_coverage   Query coverage cutoff used to align short reads (default: 50) (DIAMOND)
      --dmnd_top_pct        Keep top X% of alignments for each short read (default: 1) (DIAMOND)
      --dmnd_min_score      Minimum score for short read alignment (default: 20) (DIAMOND)
      --gencode             Genetic code used for conceptual translation (default: 11) (DIAMOND)

    For Annotation:
      --noannot             If specified, disable annotation for taxonomy or function.
                            Individual annotations can also be disabled by, e.g., setting --eggnog_db false
      --taxonomic_dmnd      Database used for taxonomic annotation 
                            (default: s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd)
      --eggnog_dmnd         One of two databases used for functional annotation with eggNOG 
                            (default: s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog_proteins.dmnd)
      --eggnog_db           One of two databases used for functional annotation with eggNOG 
                            (default: s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog.db)
      

    Batchfile:
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
include './modules/preprocess' params(
    manifest: params.manifest,
    adapter_F: params.adapter_F,
    adapter_R: params.adapter_R,
    hg_index: params.hg_index,
    min_hg_align_score: params.min_hg_align_score,
)
// Import some general tasks, such as combineReads and writeManifest
include './modules/general' params(
    savereads: params.savereads,
    output_folder: output_folder
)
include './modules/assembly' params(
    output_folder: output_folder,
    phred_offset: params.phred_offset,
    centre: params.centre,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    eggnog_db: params.eggnog_db,
    eggnog_dmnd: params.eggnog_dmnd,
    taxonomic_dmnd: params.taxonomic_dmnd,
)
include './modules/alignment' params(
    dmnd_min_identity: params.dmnd_min_identity,
    dmnd_min_coverage: params.dmnd_min_coverage,
    dmnd_top_pct: params.dmnd_top_pct,
    dmnd_min_score: params.dmnd_min_score,
    gencode: params.gencode
)


workflow {
    main:

    // Phase I: Preprocessing
    if (!params.nopreprocess) {

        // Run the entire preprocessing workflow
        preprocess_wf()

        // Combine the reads by specimen name
        combineReads(preprocess_wf.out.groupTuple())

    } else {
        // If the user specified --nopreprocess, then just 
        // read the manifest and combine by specimen
        combineReads(
            read_manifest(
                file(params.manifest)
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
            r -> [r[0], r[1]]
        }
    )

    // Make a summary of every sample and write it out to --output
    countReadsSummary(
        countReads.out.collect()
    )

    // ###################################
    // # DE NOVO ASSEMBLY AND ANNOTATION #
    // ###################################

    // Assemble samples with metaSPAdes
    metaspadesAssembly(
        combineReads.out
    )

    // Annotate those contigs with Prokka
    prokkaAnnotate(
        metaspadesAssembly.out
    )

    // Combine the gene sequences across all samples
    combineCDS(
        prokkaAnnotate.out[
            0
        ].map {
            it -> it[1]
        }.collect()
    )

    // Combine genes by amino acid identity
    clusterCDS(
        combineCDS.out
    )

    // Make a DIAMOND indexed database from those gene sequences
    makeDiamondDB(
        clusterCDS.out[0]
    )

    // Align all specimens against the DIAMOND database
    diamond(
        combineReads.out,
        makeDiamondDB.out
    )

    // Filter to the most likely single alignment per query
    famli(
        diamond.out
    )

}
