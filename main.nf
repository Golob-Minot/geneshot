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
params.adapter_F = "CTGTCTCTTATACACATCT"
params.adapter_R = "CTGTCTCTTATACACATCT"
params.hg_index_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30

// Assembly options
params.gene_fasta = false
params.assembly_folder = false
params.phred_offset = 33 // spades
params.min_identity = 90 // linclust and reference genome alignment
params.min_coverage = 50 // linclust and reference genome alignment
params.dmnd_min_identity = 80 // DIAMOND
params.dmnd_min_coverage = 50 // DIAMOND
params.dmnd_top_pct = 1 // DIAMOND
params.dmnd_min_score = 20 // DIAMOND
params.gencode = 11 //DIAMOND
params.sd_mean_cutoff = 3.0 // FAMLI
params.famli_batchsize = 10000000 // FAMLI
params.famli_folder = false // Import FAMLI-filtered alignments

// Annotation options
params.gene_shard_size = 100000
params.noannot = false
params.taxonomic_dmnd = false
params.tax_evalue = 0.00001
params.ncbi_taxdump = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
params.eggnog_db = false
params.eggnog_dmnd = false
params.eggnog_evalue = 0.00001

// CAG options
params.distance_threshold = 0.25
params.distance_metric = "cosine"
params.min_contig_size = 3
params.min_contig_depth = 5
params.min_specimens = 1
params.cags_csv = false

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
      --adapter_F           Forward sequencing adapter sequence (to be removed)
      --adapter_R           Reverse sequencing adapter sequence (to be removed)
                              (Adapter sequences default to nextera adapters)
      --min_hg_align_score  Minimum alignment score for human genome (default 30)

    For Assembly:
      --gene_fasta          (optional) Compressed FASTA with pre-generated catalog of microbial genes.
                            If provided, along with --assembly_folder then the entire de novo assembly process will be omitted.
      --assembly_folder     (optional) Folder containing de novo assemblies of every specimen.
                            If provided along with --gene_fasta, then the entire de novo assembly process will be omitted.
      --cags_csv            (optional) If a pre-generated --gene_fasta is provided, CAG groups must be provided
                            in "*.csv.gz" format from the geneshot run which generated that gene catalog.
      --phred_offset        for spades. Default 33.
      --min_identity        Amino acid identity cutoff used to combine similar genes (default: 90)
      --min_coverage        Length cutoff used to combine similar genes (default: 50) (linclust)

    For Annotation:
      --noannot             If specified, disable annotation for taxonomy or function.
                            Individual annotations can also be disabled by, e.g., setting --eggnog_db false
      --gene_shard_size     How many genes to include in a shard for annotation. (default: 100000)
      --taxonomic_dmnd      Database used for taxonomic annotation (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd)
      --tax_evalue          Maximum E-value threshold for taxonomic annotation by DIAMOND alignment
      --ncbi_taxdump        Reference describing the NCBI Taxonomy
                            (default: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
      --eggnog_dmnd         One of two databases used for functional annotation with eggNOG (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-06-17-eggNOG-v5.0/eggnog_proteins.dmnd)
      --eggnog_db           One of two databases used for functional annotation with eggNOG (default: false)
                            (Data available at s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-06-17-eggNOG-v5.0/eggnog.db)
      --eggnog_evalue       Maximum E-value threshold for eggNOG ortholog alignment (default: 0.00001)
    
    For Alignment:
      --dmnd_min_identity   Amino acid identity cutoff used to align short reads (default: 90) (DIAMOND)
      --dmnd_min_coverage   Query coverage cutoff used to align short reads (default: 50) (DIAMOND)
      --dmnd_top_pct        Keep top X% of alignments for each short read (default: 1) (DIAMOND)
      --dmnd_min_score      Minimum score for short read alignment (default: 20) (DIAMOND)
      --gencode             Genetic code used for conceptual translation (default: 11) (DIAMOND)
      --sd_mean_cutoff      Ratio of standard deviation / mean depth of sequencing used to filter genes (default: 3.0) (FAMLI)
      --famli_batchsize     Number of alignments to deduplicate in batches (default: 10000000) (FAMLI)
      --famli_folder        Optional: Specify a folder containing previously-computed FAMLI outputs

    For CAGs:
      --distance_metric     Distance metric used to group genes by co-abundance (default: cosine)
      --distance_threshold  Distance threshold used to group genes by co-abundance (default: 0.25)
      --min_contig_size     Only cluster genes which are found on contigs with at least this number of genes (default: 3; 0 to disable)
      --min_contig_depth    Minimum depth of sequencing per contig (default: 5; 0 to disable)
      --min_specimens       Only cluster genes which are detected by alignment in at least this number of specimens (default: 1)
      
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
include { preprocess_wf } from './modules/preprocess' params(
    manifest: params.manifest,
    adapter_F: params.adapter_F,
    adapter_R: params.adapter_R,
    hg_index: params.hg_index,
    hg_index_url: params.hg_index_url,
    min_hg_align_score: params.min_hg_align_score,
)
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
} from './modules/general' params(
    output_folder: output_folder,
    output_prefix: params.output_prefix,
    formula: params.formula,
    distance_metric: params.distance_metric,
    distance_threshold: params.distance_threshold,
    sd_mean_cutoff: params.sd_mean_cutoff,
    min_identity: params.min_identity,
    min_coverage: params.min_coverage,
    dmnd_min_identity: params.dmnd_min_identity,
    dmnd_min_coverage: params.dmnd_min_coverage,
    eggnog_evalue: params.eggnog_evalue,
    savereads: params.savereads,
    fdr_method: params.fdr_method,
)

// Import the workflows used for assembly
include { 
    assembly_wf;
    annotation_wf;
} from './modules/assembly' params(
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
    tax_evalue: params.tax_evalue,
    gene_fasta: params.gene_fasta,
    assembly_folder: params.assembly_folder,
)

// Import the workflows used for alignment-based analysis
include { 
    alignment_wf;
    import_alignments_wf
} from './modules/alignment' params(
    output_folder: output_folder,
    dmnd_min_identity: params.dmnd_min_identity,
    dmnd_min_coverage: params.dmnd_min_coverage,
    dmnd_top_pct: params.dmnd_top_pct,
    dmnd_min_score: params.dmnd_min_score,
    gencode: params.gencode,
    sd_mean_cutoff: params.sd_mean_cutoff,
    famli_batchsize: params.famli_batchsize,
)

// Import the workflow used to make CAGs
include { makeCAGs } from './modules/make_cags' params(
    distance_metric: params.distance_metric,
    distance_threshold: params.distance_threshold,
    min_contig_size: params.min_contig_size,
    min_contig_depth: params.min_contig_depth,
    min_specimens: params.min_specimens,
)

// Import the workflows used for statistical analysis
// Use separate workflows for corncob, each
// to analyze either the CAGs, tax IDs, or eggNOG groups
include { 
    validation_wf; 
    breakaway;
    collectBreakaway;
    corncob_wf as cag_corncob_wf;
    corncob_wf as tax_corncob_wf;
    corncob_wf as func_corncob_wf;
} from './modules/statistics' params(
    output_folder: output_folder,
    formula: params.formula,
    corncob_batches: params.corncob_batches,
    fdr_method: params.fdr_method,
    output_prefix: params.output_prefix,    
)

// Import the workflow used for composition analysis
include { metaphlan2_fastq } from './modules/composition' params(
    output_folder: output_folder
)

// Process to publish specific output files
include {
    addMetaPhlAn2Results;
    publish as publishGeneAbundances
} from './modules/general' params(
    output_folder: "${output_folder}/abund/"
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

    // ###################################
    // # DE NOVO ASSEMBLY AND ANNOTATION #
    // ###################################

    // Run the assembly and annotation workflow (in modules/assembly.nf)
    assembly_wf(
        combineReads.out,
    )

    // Run the annotation steps on the gene catalog
    annotation_wf(
        assembly_wf.out.gene_fasta
    )

    // ############################
    // # ALIGNMENT-BASED ANALYSIS #
    // ############################

    // If the user specified a folder of existing FAMLI results
    if ( params.famli_folder ) {

        // Import those results directly
        import_alignments_wf(
            params.famli_folder,
            params.output_prefix
        )

        // Point to the output of that sub-workflow
        alignments_output = import_alignments_wf.out

    // Otherwise
    } else {

        // Run the alignment-based analysis steps (in modules/alignment.nf)
        alignment_wf(
            assembly_wf.out.gene_fasta,
            combineReads.out,
            params.output_prefix
        )

        // Point to the output of the alignment sub-workflow
        alignments_output = alignment_wf.out
        
    }

    // #############
    // # MAKE CAGS #
    // #############

    // If a set of CAG assignments were provided
    if ( params.cags_csv ) {

        // Point to that file
        cags_csv = file(params.cags_csv)

    // Otherwise
    } else {

        // Make CAGs using both co-assembly and co-abundance information
        makeCAGs(
            assembly_wf.out.detailed_hdf,
            alignments_output.detailed_hdf,
        )

        // Point to the outputs from that process
        cags_csv = makeCAGs.out[0]

    }

    // ###################################
    // # STATISTICAL ANALYSIS - RICHNESS #
    // ###################################

    // Calculate the richness of each sample using the breakaway algorithm
    breakaway(
        alignments_output.famli_json_list.flatten()
    )
    collectBreakaway(
        breakaway.out.toSortedList()
    )

    // ###################
    // # GATHER RESULTS #
    // ###################

    // Start by gathering all of the results which are generated
    // no matter what options were selected by the user
    // NOTE: The code used here is imported from ./modules/general.nf

    // Join the detailed results from assembly and annotation
    joinDetailedHDF(
        assembly_wf.out.detailed_hdf,
        alignments_output.detailed_hdf,
    )

    collectResults(
        cags_csv,
        countReadsSummary.out,
        manifest_file,
        alignments_output.specimen_gene_count_csv,
        alignments_output.specimen_reads_aligned_csv,
        alignments_output.gene_length_csv,
        collectBreakaway.out,
        assembly_wf.out.n_genes_assembled_csv
    )
    resultsHDF = collectResults.out

    // ###################
    // # SPLIT CAG FASTA #
    // ###################

    splitCagFasta(
        resultsHDF,
        assembly_wf.out.gene_fasta
    )

    // ##################################
    // # STATISTICAL ANALYSIS - FORMULA #
    // ##################################

    // Calculate the association of individual CAGs with user-provided features
    cag_corncob_wf(
        resultsHDF,
        joinDetailedHDF.out,
        manifest_file,
        "cag",
        "CAG"
    )
    resultsHDF = cag_corncob_wf.out

    // ##########################
    // # COMPOSITIONAL ANALYSIS #
    // ##########################

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
                    annotation_wf.out.eggnog_tsv
                )

                // Calculate the association of functional groups with user-provided features
                func_corncob_wf(
                    addEggnogResults.out,
                    joinDetailedHDF.out,
                    manifest_file,
                    "func",
                    "eggNOG_desc_ix"
                )
                resultsHDF = func_corncob_wf.out
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
                    annotation_wf.out.tax_tsv,
                    readTaxonomy.out
                )

                // Calculate the association of functional groups with user-provided features
                tax_corncob_wf(
                    addTaxResults.out,
                    joinDetailedHDF.out,
                    manifest_file,
                    "taxa",
                    "tax_id"
                )
                resultsHDF = tax_corncob_wf.out
            }
        }
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
        joinDetailedHDF.out
    )

    // Compile a redis store with results formatted for rapid visualization
    buildRedis(
        repackFullHDF.out,
        repackDetailedHDF.out
    )

}
