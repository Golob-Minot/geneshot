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

// Annotation options
params.noannot = false
params.taxonomic_dmnd = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd"
params.ncbi_taxdump = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
params.eggnog_db = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog.db"
params.eggnog_dmnd = "s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.eggnog_proteins.dmnd"

// CAG options
params.distance_threshold = 0.25
params.distance_metric = "cosine"
params.linkage_type = "average"

// Statistical analysis options
params.formula = false

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
                            If provided, then the entire de novo assembly process will be skipped entirely.
      --phred_offset        for spades. Default 33.
      --min_identity        Amino acid identity cutoff used to combine similar genes (default: 90)
      --min_coverage        Length cutoff used to combine similar genes (default: 50) (linclust)

    For Annotation:
      --noannot             If specified, disable annotation for taxonomy or function.
                            Individual annotations can also be disabled by, e.g., setting --eggnog_db false
      --taxonomic_dmnd      Database used for taxonomic annotation 
                            (default: s3://fh-ctr-public-reference-data/tool_specific_data/geneshot/2020-01-15-geneshot/DB.refseq.tax.dmnd)
      --ncbi_taxdump        Reference describing the NCBI Taxonomy
                            (default: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
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
      --sd_mean_cutoff      Ratio of standard deviation / mean depth of sequencing used to filter genes (default: 3.0) (FAMLI)
      --famli_batchsize     Number of alignments to deduplicate in batches (default: 10000000) (FAMLI)

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
include collectAbundances from './modules/general' params(
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
include writeManifest from './modules/general' params(
    savereads: params.savereads,
    output_folder: output_folder
)
include combineReads from './modules/general' params(
    savereads: params.savereads,
    output_folder: output_folder
)
include addGeneAssembly from './modules/general'
include readTaxonomy from './modules/general'
include addEggnogResults from './modules/general'
include addCorncobResults from './modules/general'
include addTaxResults from './modules/general'
include makeSummaryHDF from './modules/general' params(
    output_prefix: params.output_prefix
)
include repackHDF as repackFullHDF from './modules/general'
include repackHDF as repackSummaryHDF from './modules/general'

// Import the workflows used for assembly
include assembly_wf from './modules/assembly' params(
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

// Import the workflows used for annotation
include annotation_wf from './modules/assembly' params(
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

<<<<<<< HEAD
// Import the workflows used for alignment-based analysis
include alignment_wf from './modules/alignment' params(
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
    famli_batchsize: params.famli_batchsize
)
=======
// Concatenate reads by sample name
process concatenate {
  container "ubuntu:16.04"
  errorStrategy "retry"
  
  input:
  set sample_name, file(fastq_list) from concatenate_ch.groupTuple()
  
  output:
  set sample_name, file("${sample_name}.fastq.gz") into correct_headers_ch

  afterScript "rm *"

  """
set -e

for fp in ${fastq_list}; do
    echo "Checking for \$fp"
    [[ -s \$fp ]]
done

cat ${fastq_list} > TEMP && mv TEMP ${sample_name}.fastq.gz
  """
>>>>>>> master

// Import the workflows used for statistical analysis
include validation_wf from './modules/statistics' params(
    output_folder: output_folder,
    formula: params.formula
)
include corncob_wf from './modules/statistics' params(
    output_folder: output_folder,
    formula: params.formula
)

<<<<<<< HEAD
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
=======
// Make sure that every read has a unique name
process correctHeaders {
  container "ubuntu:16.04"
  errorStrategy "retry"
  
  input:
  set sample_name, file(fastq) from correct_headers_ch
  
  output:
  set sample_name, file("${sample_name}.unique.headers.fastq.gz") into count_reads, metaphlan_ch, diamond_ch, humann_ch

  afterScript "rm *"

  """
set -e

gunzip -c ${fastq} | \
awk '{if(NR % 4 == 1){print("@" 1 + ((NR - 1) / 4))}else{print}}' | \
gzip -c > \
${sample_name}.unique.headers.fastq.gz
  """
>>>>>>> master

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

<<<<<<< HEAD
    }
=======
// Count the number of input reads
process countReads {
  container "ubuntu:16.04"
  errorStrategy "retry"
  
  input:
  set sample_name, file(fastq) from count_reads
  
  output:
  file "${sample_name}.countReads.csv" into total_counts
>>>>>>> master

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

    // ###################################
    // # DE NOVO ASSEMBLY AND ANNOTATION #
    // ###################################

    // A gene catalog was provided, so skip de novo assembly
    if ( params.gene_fasta ) {

<<<<<<< HEAD
        // Point to the file provided
        gene_fasta = file(params.gene_fasta)
=======
// Make a single file which summarizes the number of reads across all samples
// This is only run after all of the samples are done processing through the
// 'total_counts' channel, which is transformed by the .toSortedList() command
// into a single list containing all of the data from all samples.
process countReadsSummary {
  container "ubuntu:16.04"
  // The output from this process will be copied to the --output_folder specified by the user
  publishDir "${params.output_folder}"
  errorStrategy "retry"

  input:
  // Because the input channel has been collected into a single list, this process will only be run once
  file readcount_csv_list from total_counts.toSortedList()
  val output_prefix from params.output_prefix
  
  output:
  file "${output_prefix}.readcounts.csv" into readcounts_csv

  afterScript "rm *"

  """
set -e

echo name,n_reads > TEMP

for fp in ${readcount_csv_list}; do

    echo "Checking to make sure that \$fp exists"
    [[ -s \$fp ]]

done

cat ${readcount_csv_list} >> TEMP && mv TEMP ${output_prefix}.readcounts.csv
  """
>>>>>>> master

    } else {

<<<<<<< HEAD
        // Run the assembly and annotation workflow (in modules/assembly.nf)
        assembly_wf(
            combineReads.out
        )
=======
// Process to quantify the microbial species present using the metaphlan2 tool
process metaphlan2 {
    container "quay.io/fhcrc-microbiome/metaphlan@sha256:51b416458088e83d0bd8d840a5a74fb75066b2435d189c5e9036277d2409d7ea"
>>>>>>> master

        gene_fasta = assembly_wf.out.gene_fasta
    }

    // Run the annotation steps on the gene catalog
    annotation_wf(
        gene_fasta
    )

    // ############################
    // # ALIGNMENT-BASED ANALYSIS #
    // ############################

<<<<<<< HEAD
    // Run the alignment-based analysis steps (in modules/alignment.nf)
    alignment_wf(
        gene_fasta,
        combineReads.out
    )
=======
// If the user specifies --humann, run the tasks below
if (params.humann) {

  // Download the HUMAnN2 reference database
  process HUMAnN2_DB {
    container "quay.io/fhcrc-microbiome/humann2:v0.11.2--1"

    output:
    file "HUMANn2_DB.tar" into humann_db

    afterScript "rm -rf *"

    """
set -e

# Make a folder for the database files
mkdir HUMANn2_DB

# Download the databases
humann2_databases --download chocophlan full HUMANn2_DB
humann2_databases --download uniref uniref90_diamond HUMANn2_DB

# Tar up the database
tar cvf HUMANn2_DB.tar HUMANn2_DB

    """
  }

  // Make a channel which links the sample name to the metaphlan output
  // so that it can be combined with the raw reads for the HUMAnN2 execution
  metaphlan_for_humann
    .map{ mpn -> tuple(mpn.name.replaceFirst(/.metaphlan.tsv/, ""), mpn) }
    .set{ keyed_metaphlan_for_humann }

  // Run HUMAnN2
  process HUMAnN2 {
    container "quay.io/fhcrc-microbiome/humann2:v0.11.2--1"

    // The .join() call below combines the FASTQ and metaphlan output which share the same sample name
    input:
    set sample_name, file(fastq), file(metaphlan_output) from humann_ch.join(keyed_metaphlan_for_humann)
    val threads from 16
    // Reference database from above
    file humann_db

    output:
    set file("${sample_name}_genefamilies.tsv"), file("${sample_name}_pathabundance.tsv"), file("${sample_name}_pathcoverage.tsv") into humann_summary

    """
set -e

# Untar the database
tar xvf ${humann_db}

# Folder for output
mkdir output

humann2 \
  --input ${fastq} \
  --output output \
  --nucleotide-database HUMANn2_DB/chocophlan \
  --protein-database HUMANn2_DB/uniref \
  --threads ${threads} \
  --taxonomic-profile ${metaphlan_output}

mv output/*_genefamilies.tsv ${sample_name}_genefamilies.tsv
mv output/*_pathabundance.tsv ${sample_name}_pathabundance.tsv
mv output/*_pathcoverage.tsv ${sample_name}_pathcoverage.tsv
    """
  }

  // Summarize all of the HUMAnN2 results
  process HUMAnN2summary {
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    publishDir "${params.output_folder}"

    input:
    file humann_tsv_list from humann_summary.flatten().collect()
    val output_prefix from params.output_prefix

    output:
    file "${output_prefix}.HUMAnN2.genefamilies.csv" into humann_genefamilies_csv
    file "${output_prefix}.HUMAnN2.pathabundance.csv" into humann_pathabundance_csv
    file "${output_prefix}.HUMAnN2.pathcoverage.csv" into humann_pathcoverage_csv

    """
#!/usr/bin/env python3
import logging
import os
import pandas as pd

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [HUMAnN2summary] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

def combine_outputs(suffix, header):
    all_dat = []
    for fp in os.listdir("."):
        if fp.endswith(suffix):
            logging.info("Reading in %s" % (fp))
            d = pd.read_csv(
                fp, 
                sep="\\t", 
                comment="#",
                names=header
            )
            d["sample"] = fp.replace(suffix, "")
            all_dat.append(d)
    logging.info("Concatenating all data")
    return pd.concat(all_dat)

combine_outputs(
    "_genefamilies.tsv",
    ["gene_family", "RPK"]
).to_csv("${output_prefix}.HUMAnN2.genefamilies.csv")
logging.info("Wrote out %s" % ("${output_prefix}.HUMAnN2.genefamilies.csv"))

combine_outputs(
    "_pathabundance.tsv",
    ["pathway", "abund"]
).to_csv("${output_prefix}.HUMAnN2.pathabundance.csv")
logging.info("Wrote out %s" % ("${output_prefix}.HUMAnN2.pathabundance.csv"))

combine_outputs(
    "_pathcoverage.tsv",
    ["pathway", "cov"]
).to_csv("${output_prefix}.HUMAnN2.pathcoverage.csv")
logging.info("Wrote out %s" % ("${output_prefix}.HUMAnN2.pathcoverage.csv"))

    """
  }
>>>>>>> master

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

<<<<<<< HEAD
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
        countReadsSummary.out,
        manifest_file
    )
=======
// Align each sample against the reference database of genes using DIAMOND
process diamond {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    errorStrategy "retry"
    
    input:
    set val(sample_name), file(input_fastq) from diamond_ch
    file refdb from file(params.ref_dmnd)
    val min_id from 90
    val query_cover from 50
    val cpu from 32
    val top from 1
    val min_score from 20
    val blocks from 15
    val query_gencode from 11

    output:
    set sample_name, file("${sample_name}.aln.gz") into aln_ch

    afterScript "rm *"

    """
    set -e
    diamond \
      blastx \
      --query ${input_fastq} \
      --out ${sample_name}.aln.gz \
      --threads ${cpu} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --min-score ${min_score} \
      --query-cover ${query_cover} \
      --id ${min_id} \
      --top ${top} \
      --block-size ${blocks} \
      --query-gencode ${query_gencode} \
      --compress 1 \
      --unal 0
    """
>>>>>>> master

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

<<<<<<< HEAD
    // If we performed statistical analysis, add the results to the HDF5
    if ( params.formula ) {
        addCorncobResults(
            finalHDF,
            corncob_wf.out
        )
=======
// Filter the alignments with the FAMLI algorithm
process famli {
    container "quay.io/fhcrc-microbiome/famli@sha256:241a7db60cb735abd59f4829e8ddda0451622b6eb2321f176fd9d76297d8c9e7"
    errorStrategy "retry"
    
    input:
    set sample_name, file(input_aln) from aln_ch
    val cpu from 16
    val batchsize from 50000000

    output:
    file "${sample_name}.json.gz" into famli_json_for_summary

    afterScript "rm *"

    """
    set -e
    famli \
      filter \
      --input ${input_aln} \
      --output ${sample_name}.json \
      --threads ${cpu} \
      --batchsize ${batchsize}
    gzip ${sample_name}.json
    """
>>>>>>> master

        finalHDF = addCorncobResults.out
    }

<<<<<<< HEAD
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
=======
// Summarize all of the results from the experiment
process summarizeExperiment {
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    publishDir "${params.output_folder}"

    input:
    file metaphlan_tsv_list from metaphlan_for_summary.toSortedList()
    file famli_json_list from famli_json_for_summary.toSortedList()
    file ref_hdf5 from file(params.ref_hdf5)
    file batchfile from file(params.batchfile)
    file readcounts_csv
    val output_prefix from params.output_prefix

    output:
    file "${output_prefix}.hdf5" into output_hdf
    file "${output_prefix}.*.csv"

    afterScript "rm *"

"""
#!/usr/bin/env python3
import gzip
import json
import logging
import os
import pandas as pd

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [assembleAbundances] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Make sure that all input files are present
for fp_list in [
    "${metaphlan_tsv_list}".split(" "),
    "${famli_json_list}".split(" "),
    ["${ref_hdf5}", "${batchfile}", "${readcounts_csv}"]
]:
    for fp in fp_list:
        assert os.path.exists(fp), "Expected to find %s" % fp

# Rename the reference HDF5 to use as the output HDF5
assert os.path.exists("${ref_hdf5}")
if "${ref_hdf5}" != "${output_prefix}.hdf5":
    logging.info("Renaming ${ref_hdf5} to ${output_prefix}.hdf5")
    os.rename("${ref_hdf5}", "${output_prefix}.hdf5")
assert os.path.exists("${output_prefix}.hdf5")
# Open a connection to the output HDF5
store = pd.HDFStore("${output_prefix}.hdf5", mode="a")

# Write the batchfile to "metadata"
logging.info("Reading in %s" % ("${batchfile}"))
metadata = pd.read_csv("${batchfile}", sep=",")
logging.info("Writing metadata to HDF")
metadata.to_hdf(store, "metadata")

# Parse the list of files passed in as a space-delimited string
def parse_file_list(file_list_str, suffix):
    for fp in file_list_str.split(" "):
        assert os.path.exists(fp), "Could not find file %s" % fp
        # Yield the file path, with and without the suffix removed
        yield fp.replace(suffix, ""), fp

# Read in the KEGG KO labels
kegg_ko = pd.read_hdf(store, "/groups/KEGG_KO")

# Read in the NCBI taxid labels
taxid = pd.read_hdf(store, "/groups/NCBI_TAXID").set_index("allele")["taxid"]

# Read in all of the FAMLI results
def read_famli_json(sample_name, fp):
    logging.info("Reading in %s" % (fp))
    df = pd.DataFrame(
        json.load(gzip.open(fp, "rt"))
    )
    # Add the sample name
    df["sample"] = sample_name

    # Calculate the proportional abundance
    df["prop"] = df["depth"] / df["depth"].sum()

    # Add the taxonomic annotation
    df["taxid"] = df["id"].apply(taxid.get)

    return df

allele_abund = pd.concat([
    read_famli_json(sample_name, fp)
    for sample_name, fp in parse_file_list("${famli_json_list}", ".json.gz")
])

# Write out the FAMLI results
allele_abund.to_csv("${output_prefix}.alleles.csv", sep=",", index=None)
allele_abund.to_hdf(store, "abund/alleles", format="table", data_columns=["sample", "id"], complevel=5)

# Read in all of the MetaPhlAn2 results
def read_metaphlan(sample_name, fp):
    logging.info("Reading in %s" % (fp))
    d = pd.read_csv(
        fp, 
        sep="\\t"
    ).rename(columns=dict([
        ("Metaphlan2_Analysis", "abund")
    ]))
    # Transform into a proportion
    d["abund"] = d["abund"].apply(float) / 100

    # Add the taxonomic rank
    tax_code = dict([
        ("k", "kingdom"),
        ("p", "phylum"),
        ("c", "class"),
        ("o", "order"),
        ("f", "family"),
        ("g", "genus"),
        ("s", "species"),
        ("t", "strain"),
        ("u", "unclassified")
    ])

    d["rank"] = d["#SampleID"].apply(
        lambda s: tax_code[s.split("|")[-1][0]]
    )
>>>>>>> master

    // If we performed taxonomic analysis with DIAMOND, add the results to the HDF5
    if ( params.noannot == false ) {
        if ( params.taxonomic_dmnd ) {
            if ( !file(params.taxonomic_dmnd).isEmpty() ){
                readTaxonomy(
                    file(params.ncbi_taxdump)
                )

                addTaxResults(
                    finalHDF,
                    annotation_wf.out.tax_tsv,
                    readTaxonomy.out
                )

                finalHDF = addTaxResults.out
            }
        }
    }

    // "Repack" the HDF5, which enhances space efficiency and adds GZIP compression
    repackFullHDF(
        finalHDF
    )
<<<<<<< HEAD

    // Make a smaller summary HDF5 with a subset of the total data
    makeSummaryHDF(
        repackFullHDF.out
    )

    // "Repack" and compress that summary HDF5 as well
    repackSummaryHDF(
        makeSummaryHDF.out
    )

    publish:
        corncob_results to: "${output_folder}/stats/", enabled: params.formula
        alignment_wf.out.famli_json_list to: "${output_folder}/abund/details/"
        repackFullHDF.out to: "${output_folder}", mode: "copy", overwrite: true
        repackSummaryHDF.out to: "${output_folder}", mode: "copy", overwrite: true

}
=======
    del d["#SampleID"]

    # Add the sample name
    d["sample"] = sample_name
    return d

metaphlan_abund = pd.concat([
    read_metaphlan(sample_name, fp)
    for sample_name, fp in parse_file_list("${metaphlan_tsv_list}", ".metaphlan.tsv")
])

# Write out the MetaPhlAn2 results
metaphlan_abund.to_csv("${output_prefix}.metaphlan.csv", sep=",", index=None)
metaphlan_abund.to_hdf(store, "abund/metaphlan", format="table", data_columns=["sample", "rank", "org_name"], complevel=5)

# Summarize abundance by KEGG KO
def summarize_ko_depth(sample_name, sample_allele_abund):
    logging.info("Summarizing KO abundance for %s" % (sample_name))
    sample_allele_prop = sample_allele_abund.set_index("id")["prop"]
    sample_allele_nreads = sample_allele_abund.set_index("id")["nreads"]
    sample_ko = kegg_ko.loc[
        kegg_ko["allele"].isin(set(sample_allele_abund["id"].tolist()))
    ].copy()
    sample_ko["sample"] = sample_name
    sample_ko["prop"] = sample_ko["allele"].apply(sample_allele_prop.get)
    sample_ko["nreads"] = sample_ko["allele"].apply(sample_allele_nreads.get)
    return sample_ko.groupby(["sample", "KO"])[["prop", "nreads"]].sum().reset_index()


ko_abund = pd.concat([
    summarize_ko_depth(sample_name, sample_allele_abund)
    for sample_name, sample_allele_abund in allele_abund.groupby("sample")
])

# Write out the proportional abundance by KEGG KO
ko_abund.to_csv("${output_prefix}.KEGG_KO.csv", sep=",", index=None)
ko_abund.to_hdf(store, "abund/KEGG_KO", format="table", data_columns=["sample", "ko"], complevel=5)

# Function to summarize abundances by arbitrary groups (e.g. CAGs)
def summarize_alleles_by_group(group_key, prefix="/groups/"):
    assert group_key.startswith(prefix)
    group_name = group_key.replace(prefix, "")
    logging.info("Summarizing abundance by %s" % (group_name))

    # Read in the groupings
    group_df = pd.read_hdf(store, group_key)
    # The columns are 'allele', 'gene', and 'group'
    for k in ["allele", "gene", "group"]:
        assert k in group_df.columns.values
    group_df.set_index("allele", inplace=True)

    # Assign the group keys to the allele abundance data
    group_abund = allele_abund.copy()
    for k in ["gene", "group"]:
        group_abund[k] = group_abund["id"].apply(group_df[k].get)

    # Add up the alleles to make genes
    group_abund = group_abund.groupby(["sample", "gene", "group"])["prop"].sum().reset_index()
    # Average the genes to make groups
    group_abund = group_abund.groupby(["sample", "group"])["prop"].mean().reset_index()

    # Write out the abundance table
    group_abund.to_csv("${output_prefix}.%s.csv" % (group_name), sep=",", index=None)
    group_abund.to_hdf(store, "abund/%s" % (group_name), format="table", data_columns=["sample", "group"], complevel=5)


# Get the summary of read counts
readcounts = pd.read_csv("${readcounts_csv}")
assert "name" in readcounts.columns.values
readcounts["name"] = readcounts["name"].apply(str)
assert "n_reads" in readcounts.columns.values

# Calculate the number of aligned reads
aligned_reads = allele_abund.groupby("sample")["nreads"].sum()

# Add the column
readcounts["aligned_reads"] = readcounts["name"].apply(aligned_reads.get)
assert readcounts["aligned_reads"].isnull().sum() == 0, (readcounts.loc[readcounts["aligned_reads"].isnull()])

# Write to HDF and CSV
readcounts.to_hdf(store, "readcounts", format="table")
readcounts.to_csv("${output_prefix}.readcounts.csv", sep=",", index=None)

for key in store.keys():
    if key.startswith("/groups/") and key != "/groups/KEGG_KO" and key != "/groups/NCBI_TAXID":
        summarize_alleles_by_group(key)

"""

}

// Add the HUMAnN2 results to the summary of the experiment, if specified
if (params.humann) {
  process addHUMAnN2toHDF {
      container "quay.io/fhcrc-microbiome/python-pandas:latest"
      publishDir "${params.output_folder}"

      input:
      file humann_genefamilies_csv
      file humann_pathabundance_csv
      file humann_pathcoverage_csv
      file output_hdf

      output:
      file "${output_hdf}"

      afterScript "rm *"

"""
#!/usr/bin/env python3
import gzip
import json
import logging
import os
import pandas as pd

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [addHUMAnN2toHDF] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Open a connection to the output HDF5
store = pd.HDFStore("${output_hdf}", mode="a")

for fp, key, dtype_dict in [
    (
      "${humann_genefamilies_csv}", 
      "/abund/humann_genefamilies", 
      dict([("gene_family", str), ("RPK", float), ("sample", str)]),
    ),
    (
      "${humann_pathabundance_csv}", 
      "/abund/humann_pathabundance", 
      dict([("pathway", str), ("abund", float), ("sample", str)]),
    ),
    (
      "${humann_pathcoverage_csv}", 
      "/abund/humann_pathcoverage", 
      dict([("pathway", str), ("cov", float), ("sample", str)],
    ))
]:
    print("Reading in %s" % fp)
    df = pd.read_csv(fp, sep=",", dtype=dtype_dict, usecols=list(dtype_dict.keys()))
    print(df.head())
    df.to_hdf(store, key, complevel=5, format="table")

store.close()
"""
  }
}
>>>>>>> master
