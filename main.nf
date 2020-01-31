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
params.gene_fasta = false
params.phred_offset = 33 // spades
params.min_identity = 90 // mmseqs2 and reference genome alignment
params.min_coverage = 50 // mmseqs2 and reference genome alignment
params.dmnd_min_identity = 80 // DIAMOND
params.dmnd_min_coverage = 50 // DIAMOND
params.dmnd_top_pct = 1 // DIAMOND
params.dmnd_min_score = 20 // DIAMOND
params.gencode = 11 //DIAMOND

// Annotation options
params.noannot = false
params.taxonomic_dmnd = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd"
params.eggnog_db = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog.db"
params.eggnog_dmnd = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog_proteins.dmnd"

// CAG options
params.distance_threshold = 0.5
params.distance_metric = "cosine"
params.linkage_type = "average"

// Statistical analysis options
params.formula = false

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
      --gene_fasta          (optional) Compressed FASTA with pre-generated catalog of microbial genes.
                            If provided, then the entire de novo assembly process will be skipped entirely.
      --phred_offset        for spades. Default 33.
      --min_identity        Amino acid identity cutoff used to combine similar genes (default: 90)
                            This value is also used to align the gene catalog against whole microbial genomes.
      --min_coverage        Length cutoff used to combine similar genes (default: 50) (mmseqs2)
                            This value is also used to align the gene catalog against whole microbial genomes.

    For Annotation:
      --noannot             If specified, disable annotation for taxonomy or function.
                            Individual annotations can also be disabled by, e.g., setting --eggnog_db false
      --taxonomic_dmnd      Database used for taxonomic annotation 
                            (default: s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd)
      --eggnog_dmnd         One of two databases used for functional annotation with eggNOG 
                            (default: s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog_proteins.dmnd)
      --eggnog_db           One of two databases used for functional annotation with eggNOG 
                            (default: s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog.db)
    
    For Alignment:
      --dmnd_min_identity   Amino acid identity cutoff used to align short reads (default: 90) (DIAMOND)
      --dmnd_min_coverage   Query coverage cutoff used to align short reads (default: 50) (DIAMOND)
      --dmnd_top_pct        Keep top X% of alignments for each short read (default: 1) (DIAMOND)
      --dmnd_min_score      Minimum score for short read alignment (default: 20) (DIAMOND)
      --gencode             Genetic code used for conceptual translation (default: 11) (DIAMOND)

    For CAGs:
      --distance_metric     Distance metric used to group genes by co-abundance (default: cosine)
      --distance_threshold  Distance threshold used to group genes by co-abundance (default: 0.5)
      --linkage_type        Linkage type used to group genes by co-abundance (default: average)
      

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
    hg_index_url: params.hg_index_url,
    min_hg_align_score: params.min_hg_align_score,
)
// Import some general tasks, such as combineReads and writeManifest
include './modules/general' params(
    savereads: params.savereads,
    output_folder: output_folder
)

// Import the workflows used for assembly and annotation
include './modules/assembly' params(
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
include './modules/alignment' params(
    output_folder: output_folder,
    dmnd_min_identity: params.dmnd_min_identity,
    dmnd_min_coverage: params.dmnd_min_coverage,
    dmnd_top_pct: params.dmnd_top_pct,
    dmnd_min_score: params.dmnd_min_score,
    gencode: params.gencode,
    distance_metric: params.distance_metric,
    distance_threshold: params.distance_threshold,
    linkage_type: params.linkage_type,
)

include './modules/statistics' params(
    output_folder: output_folder,
    formula: params.formula
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

    // A gene catalog was provided, so skip de novo assembly
    if ( params.gene_fasta ) {

        // Point to the file provided
        gene_fasta = file(params.gene_fasta)

    } else {

        // Run the assembly and annotation workflow (in modules/assembly.nf)
        assembly_wf(
            combineReads.out
        )

        gene_fasta = assembly_wf.out.gene_fasta
    }

    // Run the annotation steps on the gene catalog
    annotation_wf(
        gene_fasta
    )

    // ############################
    // # ALIGNMENT-BASED ANALYSIS #
    // ############################

    // Run the alignment-based analysis steps (in modules/alignment.nf)
    alignment_wf(
        gene_fasta,
        combineReads.out
    )

    // ########################
    // # STATISTICAL ANALYSIS #
    // ########################

    // Calculate the association of individual CAGs with user-provided features
    if ( params.formula ) {
        corncob_wf(
            alignment_wf.out.famli_json_list,
            alignment_wf.out.cag_csv,
            file(params.manifest)
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
        alignment_wf.out.cag_csv,
        alignment_wf.out.gene_abund_feather,
        alignment_wf.out.cag_abund_feather,
        alignment_wf.out.famli_json_list,
        countReadsSummary.out
    )

    // If we performed de novo assembly, add the gene assembly information
    if ( params.gene_fasta ) {
        finalHDF = collectAbundances.out
    } else {
        addGeneAssembly(
            collectAbundances.out,
            assembly_wf.out.allele_gene_tsv,
            assembly_wf.out.allele_assembly_csv
        )
        finalHDF = addGeneAssembly.out
    }

    // If we performed statistical analysis, add the results to the HDF5
    if ( params.formula ) {
        addCorncobResults(
            finalHDF,
            corncob_wf.out
        )

        finalHDF = addCorncobResults.out
    }

    // If we performed functional analysis with eggNOG, add the results to the HDF5
    if ( params.noannot == false ) {
        if ( params.eggnog_db && params.eggnog_dmnd ) {
            if ( !file(params.eggnog_db).isEmpty() && !file(params.eggnog_dmnd).isEmpty() ){
                addEggnogResults(
                    finalHDF,
                    annotation_wf.out.eggnog_tsv
                )

                finalHDF = addEggnogResults.out
            }
        }
    }

    // If we performed taxonomic analysis with DIAMOND, add the results to the HDF5
    if ( params.noannot == false ) {
        if ( params.taxonomic_dmnd ) {
            if ( !file(params.taxonomic_dmnd).isEmpty() ){
                addTaxResults(
                    finalHDF,
                    annotation_wf.out.tax_tsv
                )

                finalHDF = addTaxResults.out
            }
        }
    }

    // "Repack" the HDF5, which enhances space efficiency and adds GZIP compression
    repackHDF(
        finalHDF
    )

    publish:
        corncob_results to: "${output_folder}/stats/", enabled: params.formula
        alignment_wf.out.famli_json_list to: "${output_folder}/abund/details/"
        repackHDF.out to: "${output_folder}", mode: "copy", overwrite: true

}
