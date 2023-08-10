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
params.hg_index_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30

// Assembly options
params.gene_fasta = false
params.phred_offset = 33 // spades
params.min_identity = 90 // linclust and reference genome alignment
params.min_coverage = 50 // linclust and reference genome alignment

// Alignment options
params.dmnd_min_identity = 80 // DIAMOND
params.dmnd_min_coverage = 50 // DIAMOND
params.dmnd_top_pct = 1 // DIAMOND
params.dmnd_min_score = 20 // DIAMOND
params.gencode = 11 //DIAMOND
params.sd_mean_cutoff = 3.0 // FAMLI
params.famli_batchsize = 10000000 // FAMLI

// Annotation options
params.noannot = false
params.taxonomic_dmnd = false
params.ncbi_taxdump = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
params.eggnog_db = false
params.eggnog_dmnd = false

// CAG options
params.distance_threshold = 0.25
params.distance_metric = "cosine"
params.linkage_type = "average"
params.cag_batchsize = 10000

// Statistical analysis options
params.formula = false
params.fdr_method = "fdr_bh"
params.corncob_batches = 10

// Compositional analysis options
params.composition = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run Golob-Minot/geneshot <ARGUMENTS>
    
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
      --min_hg_align_score  Minimum alignment score for human genome (default 30)

    For Assembly:
      --gene_fasta          (optional) Compressed FASTA with pre-generated catalog of microbial genes.
                            If provided, then the entire de novo assembly process will be skipped entirely.
      --phred_offset        for spades. Default 33.
      --min_identity        Amino acid identity cutoff used to combine similar genes (default: 90)
      --min_coverage        Length cutoff used to combine similar genes (default: 50) (linclust)

    For Annotation:
      --noannot             If specified, disable annotation for taxonomy or function.
                            Individual annotations can also be disabled by, e.g., setting --eggnog_db false
      --taxonomic_dmnd      Database used for taxonomic annotation (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd)
      --ncbi_taxdump        Reference describing the NCBI Taxonomy
                            (default: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
      --eggnog_dmnd         One of two databases used for functional annotation with eggNOG (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-06-17-eggNOG-v5.0/eggnog_proteins.dmnd)
      --eggnog_db           One of two databases used for functional annotation with eggNOG (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-06-17-eggNOG-v5.0/eggnog.db)
    
    For Alignment:
      --dmnd_min_identity   Amino acid identity cutoff used to align short reads (default: 90) (DIAMOND)
      --dmnd_min_coverage   Query coverage cutoff used to align short reads (default: 50) (DIAMOND)
      --dmnd_top_pct        Keep top X% of alignments for each short read (default: 1) (DIAMOND)
      --dmnd_min_score      Minimum score for short read alignment (default: 20) (DIAMOND)
      --gencode             Genetic code used for conceptual translation (default: 11) (DIAMOND)
      --sd_mean_cutoff      Ratio of standard deviation / mean depth of sequencing used to filter genes (default: 3.0) (FAMLI)
      --famli_batchsize     Number of alignments to deduplicate in batches (default: 10000000) (FAMLI)

    For CAGs:
      --distance_metric     Distance metric used to group genes by co-abundance (default: cosine)
      --distance_threshold  Distance threshold used to group genes by co-abundance (default: 0.25)
      --linkage_type        Linkage type used to group genes by co-abundance (default: average)
      
    For Statistical Analysis:
      --formula             Optional formula used to estimate associations with CAG relative abundance
      --fdr_method          FDR method used to calculate q-values for associations (default: 'fdr_bh')
      --corncob_batches     Number of parallel processes to use processing each formula

    For Compositional Analysis:
      --composition         When included, metaPhlAn2 will be run on all specimens

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
include { Preprocess_wf } from './modules/preprocess' params(
    manifest: params.manifest,
    hg_index: params.hg_index,
    hg_index_url: params.hg_index_url,
    min_hg_align_score: params.min_hg_align_score,
    savereads: params.savereads,
    output: output_folder
)
include { CombineReads} from './modules/preprocess' params(
    savereads: params.savereads,
    output: output_folder
)
include { WriteManifest} from './modules/preprocess' params(
    savereads: params.savereads,
    output: output_folder
)
// Import some general tasks, such as CombineReads and writeManifest
include { Read_manifest } from './modules/general'
include { countReads } from './modules/general'
include { countReadsSummary } from './modules/general' params(
    output_folder: output_folder
)
include { collectAbundances } from './modules/general' params(
    output_prefix: params.output_prefix,
    formula: params.formula,
    distance_metric: params.distance_metric,
    distance_threshold: params.distance_threshold,
    linkage_type: params.linkage_type,
    sd_mean_cutoff: params.sd_mean_cutoff,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    dmnd_min_identity: params.dmnd_min_identity,
    dmnd_min_coverage: params.dmnd_min_coverage
)

include { addGeneAssembly } from './modules/general'
include { readTaxonomy } from './modules/general'
include { addEggnogResults } from './modules/general'
include { addCorncobResults } from './modules/general' params(
    fdr_method: params.fdr_method
)
include { addTaxResults } from './modules/general'
include { repackHDF as repackFullHDF } from './modules/general' params(
    output_folder: output_folder
)
include { repackHDF as repackDetailedHDF } from './modules/general' params(
    output_folder: output_folder
)

// Import the workflows used for assembly
include { Genecatalog_wf } from './modules/genecatalog' params(
    output_folder: output_folder,
    output_prefix: params.output_prefix,
    phred_offset: params.phred_offset,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    noannot: params.noannot,
    eggnog_db: params.eggnog_db,
    eggnog_dmnd: params.eggnog_dmnd,
    taxonomic_dmnd: params.taxonomic_dmnd,
    gencode: params.gencode,
)

// Import the workflows used for annotation
include { Annotation_wf } from './modules/annotation' params(
    output_folder: output_folder,
    phred_offset: params.phred_offset,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    noannot: params.noannot,
    eggnog_db: params.eggnog_db,
    eggnog_dmnd: params.eggnog_dmnd,
    taxonomic_dmnd: params.taxonomic_dmnd,
    gencode: params.gencode,
)

// Import the workflows used for alignment-based analysis
include { Alignment_wf } from './modules/quantify' params(
    output_folder: output_folder,
    dmnd_min_identity: params.dmnd_min_identity,
    dmnd_min_coverage: params.dmnd_min_coverage,
    dmnd_top_pct: params.dmnd_top_pct,
    dmnd_min_score: params.dmnd_min_score,
    gencode: params.gencode,
    sd_mean_cutoff: params.sd_mean_cutoff,
    famli_batchsize: params.famli_batchsize,
    cag_batchsize: params.cag_batchsize
)

// And for CAG generation
include { CAG_contig_oriented_wf } from './modules/make_cags' params(
    output_folder: output_folder,
    dmnd_min_identity: params.dmnd_min_identity,
    dmnd_min_coverage: params.dmnd_min_coverage,
    dmnd_top_pct: params.dmnd_top_pct,
    dmnd_min_score: params.dmnd_min_score,
    gencode: params.gencode,
    distance_metric: params.distance_metric,
    distance_threshold: params.distance_threshold,
    linkage_type: params.linkage_type,
    sd_mean_cutoff: params.sd_mean_cutoff,
    famli_batchsize: params.famli_batchsize,
    cag_batchsize: params.cag_batchsize
)

// Import the workflows used for statistical analysis
include { validation_wf } from './modules/statistics' params(
    output_folder: output_folder,
    formula: params.formula,
    corncob_batches: params.corncob_batches
)
include { corncob_wf } from './modules/statistics' params(
    output_folder: output_folder,
    formula: params.formula,
    corncob_batches: params.corncob_batches
)
include { runBetta } from './modules/statistics'
include { addBetta } from './modules/statistics' params(
    fdr_method: params.fdr_method
)
include { breakaway } from './modules/statistics'
include { collectBreakaway } from './modules/statistics' params(
    output_folder: output_folder,
    output_prefix: params.output_prefix
)

// Import the workflow used for composition analysis
include { metaphlan2_fastq } from './modules/composition' params(
    output_folder: output_folder
)
// include join_metaphlan2 from './modules/composition'
include { addMetaPhlAn2Results } from './modules/general'

// Process to publish specific output files
include { publish as publishGeneAbundances } from './modules/general' params(
    output_folder: "${output_folder}/abund/"
)

workflow {
    main:

    // Phase 0: Validation of input data

    // If the user specifies a `--formula`, the first step in the process
    // will be to ensure that the formula is written correctly, and is
    // compatible with the data provided in the manifest
    if ( params.formula ) {
        // Set up a channel with the strings of the formula(s) provided
        formula_ch = Channel.of(
            params.formula.split(",")
        )
        validation_wf(
            file(params.manifest),
            formula_ch
        )
        manifest_file = validation_wf.out
    } else {
        manifest_file = Channel.from(file(params.manifest))
    }

    manifest_qced = Read_manifest(manifest_file)

    // Phase I: Preprocessing
    if (!params.nopreprocess) {

        // Run the entire preprocessing workflow
        Preprocess_wf(
            manifest_qced.valid_paired_indexed,
             manifest_qced.valid_paired
        )

        combined_reads = Preprocess_wf

    } else {
        // If the user specified --nopreprocess, then just 
        // read the manifest and combine by specimen
        CombineReads(
            manifest_qced.valid_paired.mix(manifest_qced.valid_paired_indexed)
            .map { 
                r -> [r.specimen, file(r.R1), file(r.R2)]
            }.groupTuple()
        )
        combined_reads = CombineReads
        // If the user specified --savereads, write out the manifest
        if (params.savereads) {
            writeManifest(
                combined_reads.out
            )
        }        
    }

    // Count the reads for every sample individually (just take the first of the pair of reads)
    countReads(
        combined_reads.out.map {
            r -> [r[0], r[1], r[2]]
        }
    )

    // Make a summary of every sample and write it out to --output
    countReadsSummary(
        countReads.collect()
    )

    // ##########################
    // # COMPOSITIONAL ANALYSIS #
    // ##########################
    if (params.composition) {
        metaphlan2_fastq(
            combined_reads.out.map {
                r -> [r[0], r[1], r[2]]
            }
        )
    }

    // ###################################
    // # DE NOVO ASSEMBLY AND ANNOTATION #
    // ###################################

    // A gene catalog was provided, so skip de novo assembly
    if ( params.gene_fasta ) {

        // Point to the file provided
        gene_fasta = file(params.gene_fasta)

    } else {

        // Run the assembly and annotation workflow (in modules/genecatalog.nf)
        Genecatalog_wf(
            combined_reads.out
        )

        gene_fasta = Genecatalog_wf.out.gene_fasta
    }

    // Run the annotation steps on the gene catalog
    Annotation_wf(
        gene_fasta
    )

    // ############################
    // # ALIGNMENT-BASED ANALYSIS #
    // ############################

    // Run the alignment-based analysis steps (in modules/alignment.nf)
    Alignment_wf(
        gene_fasta,
        combined_reads.out,
        'gene' // output prefix here, of gene level output
    )

    /*
    // And the CAG generation steps
    CAG_contig_oriented_wf(
        Alignment_wf.out.gene_abundances_zarr_tar,
        Alignment_wf.out.gene_lists
    )

    // ########################
    // # STATISTICAL ANALYSIS #
    // ########################

    // Calculate the richness of each sample using the breakaway algorithm
    breakaway(
        Alignment_wf.out.famli_json_list.flatten()
    )
    collectBreakaway(
        breakaway.out.toSortedList()
    )

    // Calculate the association of individual CAGs with user-provided features
    if ( params.formula ) {
        corncob_wf(
            Alignment_wf.out.famli_json_list,
            CAG_contig_oriented_wf.out.cag_csv,
            file(params.manifest),
            formula_ch
        )
        corncob_results = corncob_wf.out
    } else {
        corncob_results = Channel.empty()
    }

    // ###################
    // # GATHER RESULTS #
    // ###################

    // Start by gathering all of the results which are generated
    // no matter what options were selected by the user
    // NOTE: The code used here is imported from ./modules/general.nf

    collectAbundances(
        CAG_contig_oriented_wf.out.cag_csv,
        CAG_contig_oriented_wf.out.cag_abund_feather,
        countReadsSummary.out,
        manifest_file,
        Alignment_wf.out.specimen_gene_count_csv,
        Alignment_wf.out.specimen_reads_aligned_csv,
        Alignment_wf.out.gene_length_csv,
        collectBreakaway.out,
    )

    // If we performed de novo assembly, add the gene assembly information
    if ( params.gene_fasta ) {
        resultsHDF = collectAbundances.out
        detailedHDF = Alignment_wf.out.detailed_hdf
    } else {
        addGeneAssembly(
            collectAbundances.out,
            Alignment_wf.out.detailed_hdf,
            Genecatalog_wf.out.allele_assembly_csv_list
        )
        resultsHDF = addGeneAssembly.out[0]
        detailedHDF = addGeneAssembly.out[1]
    }

    // If we performed compositional analysis, add the results ot the HDF5
    if (params.composition) {
        addMetaPhlAn2Results(
            resultsHDF,
            metaphlan2_fastq.out.map {
                r -> r[1]
            }.toSortedList()
        )

        resultsHDF = addMetaPhlAn2Results.out
    }

    // If we performed functional analysis with eggNOG, add the results to the HDF5
    if ( params.noannot == false ) {
        if ( params.eggnog_db && params.eggnog_dmnd ) {
            if ( !file(params.eggnog_db).isEmpty() && !file(params.eggnog_dmnd).isEmpty() ){
                addEggnogResults(
                    resultsHDF,
                    Annotation_wf.out.eggnog_tsv
                )

                resultsHDF = addEggnogResults.out
            }
        }
    }

    // If we performed taxonomic analysis with DIAMOND, add the results to the HDF5
    if ( params.noannot == false ) {
        if ( params.taxonomic_dmnd ) {
            if ( !file(params.taxonomic_dmnd).isEmpty() ){
                readTaxonomy(
                    file(params.ncbi_taxdump)
                )

                addTaxResults(
                    resultsHDF,
                    Annotation_wf.out.tax_tsv,
                    readTaxonomy.out
                )

                resultsHDF = addTaxResults.out
            }
        }
    }

    // If we performed statistical analysis, add the results to the HDF5
    if ( params.formula ) {
        addCorncobResults(
            resultsHDF,
            corncob_wf.out
        )

        runBetta(
            addCorncobResults.out[1].flatten()
        )

        addBetta(
            addCorncobResults.out[0],
            runBetta.out.toSortedList()
        )

        resultsHDF = addBetta.out[0]

    }

    // "Repack" the HDF5, which enhances space efficiency and adds GZIP compression
    repackFullHDF(
        resultsHDF
    )

    // "Repack" and compress the detailed results HDF5 as well
    repackDetailedHDF(
        detailedHDF
    )
    // */

}
