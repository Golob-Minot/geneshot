#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

  The quantify_cags.nf workflow will align a set of short sequence reads against a previously-generated
  gene catalog and output the relative abundances of genes and CAGs in those samples.

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
params.output_prefix = 'geneshot'
params.manifest = null

// Preprocessing options
params.adapter_F = "CTGTCTCTTATACACATCT"
params.adapter_R = "CTGTCTCTTATACACATCT"
params.hg_index_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30

// Genes and CAGs from the previous analysis
params.gene_fasta = false
params.cags_csv = false

// Alignment parameters
params.dmnd_min_identity = 80 // DIAMOND
params.dmnd_min_coverage = 50 // DIAMOND
params.dmnd_top_pct = 1 // DIAMOND
params.dmnd_min_score = 20 // DIAMOND
params.gencode = 11 //DIAMOND
params.sd_mean_cutoff = 3.0 // FAMLI
params.famli_batchsize = 10000000 // FAMLI
params.famli_folder = false // Import FAMLI-filtered alignments

// Compositional analysis options
params.composition = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run Golob-Minot/geneshot/quantify_cags.nf <ARGUMENTS>
    
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

    Gene catalog:
      --gene_fasta          (required) Compressed FASTA with pre-generated catalog of microbial genes.
      --cags_csv            (required) CAG groups must be provided in "*.csv.gz" format from the geneshot run which generated that gene catalog.

    For Alignment:
      --dmnd_min_identity   Amino acid identity cutoff used to align short reads (default: 90) (DIAMOND)
      --dmnd_min_coverage   Query coverage cutoff used to align short reads (default: 50) (DIAMOND)
      --dmnd_top_pct        Keep top X% of alignments for each short read (default: 1) (DIAMOND)
      --dmnd_min_score      Minimum score for short read alignment (default: 20) (DIAMOND)
      --gencode             Genetic code used for conceptual translation (default: 11) (DIAMOND)
      --sd_mean_cutoff      Ratio of standard deviation / mean depth of sequencing used to filter genes (default: 3.0) (FAMLI)
      --famli_batchsize     Number of alignments to deduplicate in batches (default: 10000000) (FAMLI)
      --famli_folder        Optional: Specify a folder containing previously-computed FAMLI outputs

    For Compositional Analysis:
      --composition         When included, metaPhlAn2 will be run on all specimens

    Manifest:
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
include { preprocess_wf } from './modules/preprocess'

// Import some general tasks, such as combineReads and writeManifest
include { 
    countReads;
    read_manifest;
    countReadsSummary;
    collectResults;
    writeManifest;
    combineReads;
    readTaxonomy;
    addEggnogResults;
    addTaxResults;
    calcDistances;
    joinHDF as joinDetailedHDF;
    repackHDF as repackFullHDF;
    repackHDF as repackDetailedHDF;
    buildRedis;
    splitCagFasta
} from './modules/general' addParams(
    output_folder: output_folder
)

// Import the workflows used for alignment-based analysis
include { 
    alignment_wf;
    import_alignments_wf
} from './modules/alignment' addParams(
    output_folder: output_folder
)

// Import the workflow used to make CAGs
include { makeCAGs } from './modules/make_cags'

// Import the workflows used for statistical analysis
// Use separate workflows for corncob, each
// to analyze either the CAGs, tax IDs, or eggNOG groups
include { 
    breakaway;
    collectBreakaway
} from './modules/statistics' addParams(
    output_folder: output_folder
)

// Import the workflow used for composition analysis
include { metaphlan2_fastq } from './modules/composition' addParams(
    output_folder: output_folder
)

// Process to publish specific output files
include {
    addMetaPhlAn2Results;
    publish as publishGeneAbundances
} from './modules/general' addParams(
    output_folder: "${output_folder}/abund/"
)

workflow {
    main:

    log.info"""

    Running: Golob-Minot/geneshot/quantify_cags.nf
    Working files: ${workDir}
    
    Required parameters:
      manifest:            "${params.manifest}"
      gene_fasta:          "${params.gene_fasta}"
      cags_csv:            "${params.cags_csv}"
      output:              "${params.output}"
      output_prefix:       "${params.output_prefix}"

    Optional parameters
      nopreprocess:        "${params.nopreprocess}"
      savereads:           "${params.savereads}"
      hg_index_url:        "${params.hg_index_url}"
      hg_index:            "${params.hg_index}"
      adapter_F:           "${params.adapter_F}"
      adapter_R:           "${params.adapter_R}"
      min_hg_align_score:  "${params.min_hg_align_score}"
      dmnd_min_identity:   "${params.dmnd_min_identity}"
      dmnd_min_coverage:   "${params.dmnd_min_coverage}"
      dmnd_top_pct:        "${params.dmnd_top_pct}"
      dmnd_min_score:      "${params.dmnd_min_score}"
      gencode:             "${params.gencode}"
      sd_mean_cutoff:      "${params.sd_mean_cutoff}"
      famli_batchsize:     "${params.famli_batchsize}"
      famli_folder:        "${params.famli_folder}"
      composition:         "${params.composition}"
    """

    if(!params.gene_fasta){
        error "Must provide --gene_fasta"
    }

    if(!params.cags_csv){
        error "Must provide --cags_csv"
    }

    if(!params.output){
        error "Must provide --output"
    }

    // Phase 0: Validation of input data

    // Read the input data
    manifest_file = Channel.from(file(params.manifest, checkIfExists: true))
    manifest_qced = read_manifest(manifest_file)

    // Phase I: Preprocessing
    if (!params.nopreprocess) {

        // Run the entire preprocessing workflow
        preprocess_wf(
            manifest_qced.valid_paired_indexed,
            manifest_qced.valid_paired
        )

        // Combine the reads by specimen name
        combineReads(preprocess_wf.out.groupTuple())

    } else {
        // If the user specified --nopreprocess, then just 
        // read the manifest and combine by specimen
        combineReads(
            manifest_qced.valid_paired.mix(manifest_qced.valid_paired_indexed)
            .map { 
                r -> [r.specimen, file(r.R1), file(r.R2)]
            }.groupTuple()
        )
    }

    // If no reads were found, raise an error
    combineReads.out.ifEmpty { error "No reads found" }

    // If the user specified --savereads, write out the manifest
    if (params.savereads) {
        writeManifest(
            combineReads.out
        )
    }

    // Count the reads for every sample individually (just take the first of the pair of reads)
    countReads(
        combineReads
            .out
            .ifEmpty{ error "No reads found" }
            .map {
                r -> [r[0], r[1], r[2]]
            }
    )

    // Make a summary of every sample and write it out to --output
    countReadsSummary(
        countReads.out.collect()
    )

    // ##########################
    // # COMPOSITIONAL ANALYSIS #
    // ##########################
    if (params.composition) {
        metaphlan2_fastq(
            combineReads.out.map {
                r -> [r[0], r[1], r[2]]
            }
        )
    }

    // Run the alignment-based analysis steps (in modules/alignment.nf)
    gene_fasta = file("${params.gene_fasta}", checkIfExists: true)
    alignment_wf(
        gene_fasta,
        combineReads.out,
        params.output_prefix
    )

    // Point to the output of the alignment sub-workflow
    alignments_output = alignment_wf.out

    // Point to the CAGs CSV file
    cags_csv = file("${params.cags_csv}", checkIfExists: true)

    // ###################################
    // # STATISTICAL ANALYSIS - RICHNESS #
    // ###################################

    // Calculate the richness of each sample using the breakaway algorithm
    breakaway(
        alignments_output
            .famli_json_list
            .flatten()
    )
    collectBreakaway(
        breakaway
            .out
            .ifEmpty{ error "No alignments found" }
            .toSortedList()
    )

    // ###################
    // # GATHER RESULTS #
    // ###################

    // Start by gathering all of the results which are generated
    // no matter what options were selected by the user
    // NOTE: The code used here is imported from ./modules/general.nf

    collectResults(
        cags_csv,
        countReadsSummary.out,
        manifest_file,
        alignments_output.specimen_gene_count_csv,
        alignments_output.specimen_reads_aligned_csv,
        alignments_output.gene_length_csv,
        collectBreakaway.out,
        file("$projectDir/templates/gene_count_template.csv")
    )
    resultsHDF = collectResults.out

    // If we performed compositional analysis, add the results to the HDF5
    if (params.composition) {
        addMetaPhlAn2Results(
            resultsHDF,
            metaphlan2_fastq.out.map {
                r -> r[1]
            }.toSortedList()
        )

        resultsHDF = addMetaPhlAn2Results.out
    }

    // Add the pairwise distances between every specimen
    calcDistances(
        resultsHDF
    )

    // "Repack" the HDF5, which enhances space efficiency and adds GZIP compression
    repackFullHDF(
        calcDistances.out
    )

    // "Repack" and compress the detailed results HDF5 as well
    repackDetailedHDF(
        alignments_output.detailed_hdf
    )

    // Compile a redis store with results formatted for rapid visualization
    buildRedis(
        repackFullHDF.out,
        repackDetailedHDF.out
    )

}
