#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// containers
container__FAMLI = 'quay.io/fhcrc-microbiome/famli:v1.5'
container__anndata = 'golob/python-anndata:0.9.2'

// Processes used for alignment of reads against gene databases


// Alignment options
params.dmnd_min_identity = 80 // DIAMOND
params.dmnd_min_coverage = 50 // DIAMOND
params.dmnd_top_pct = 1 // DIAMOND
params.dmnd_min_score = 20 // DIAMOND
params.gencode = 11 //DIAMOND
params.sd_mean_cutoff = 3.0 // FAMLI
params.famli_batchsize = 10000000 // FAMLI



workflow Alignment_wf {
    take:
        gene_fasta
        reads_ch
        output_prefix

    main:

    // Make a DIAMOND indexed database from those gene sequences
    DiamondDB(
        gene_fasta
    )

    // Align all specimens against the DIAMOND database
    Diamond(
        reads_ch,
        DiamondDB.out
    )

    // Filter to the most likely single alignment per query
    Famli(
        Diamond.out
    )

    AssembleAnndata(
        Famli.out.toSortedList(),
        gene_fasta,
    )

    /*
    // Make a single table with the abundance of every gene across every sample
    AssembleAbundances(
        Famli.out.toSortedList(),
        params.cag_batchsize,
        output_prefix
    )

    // Output the abundances
    OutputSpecimenGeneAbundance(
        AssembleAbundances.out.abundance_zarr_tar
    )
    emit:
        famli_json_list = Famli.out.toSortedList()
        specimen_gene_count_csv = AssembleAbundances.out.specimen_gene_count_csv
        specimen_reads_aligned_csv = AssembleAbundances.out.specimen_reads_aligned_csv
        detailed_hdf = AssembleAbundances.out.detailed_hdf
        gene_length_csv = AssembleAbundances.out.gene_length_csv
        gene_abundances_zarr_tar = AssembleAbundances.out.abundance_zarr_tar
        gene_lists = AssembleAbundances.out.gene_lists
    // */
}


// Align each sample against the reference database of genes using DIAMOND
process DiamondDB {
    tag "Make a DIAMOND database"
    container "quay.io/fhcrc-microbiome/famli:v1.5"
    label 'mem_veryhigh'
    errorStrategy 'finish'
    publishDir "${params.output_folder}/ref/", mode: "copy"
    
    input:
    file fasta

    output:
    file "genes.dmnd"

    """
    set -e
    diamond \
      makedb \
      --in ${fasta} \
      --db genes.dmnd \
      --threads ${task.cpus}
    """

}


// Align each sample against the reference database of genes using DIAMOND
process Diamond {
    tag "Align to the gene catalog"
    container "${container__FAMLI}"
    label 'mem_veryhigh'
    errorStrategy 'finish'
    
    input:
    tuple val(sample_name), file(R1), file(R2)
    file refdb
    
    output:
    tuple val(sample_name), file("${sample_name}.aln.gz")

    """
    set -e

    for fp in ${R1} ${R2}; do
        cat \$fp | \
        gunzip -c | \
        awk "{if(NR % 4 == 1){print \\"@\$fp\\" NR }else{if(NR % 4 == 3){print \\"+\\"}else{print}}}" | \
        gzip -c \
        > query.fastq.gz

        diamond \
        blastx \
        --query query.fastq.gz \
        --out \$fp.aln.gz \
        --threads ${task.cpus} \
        --db ${refdb} \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
        --min-score ${params.dmnd_min_score} \
        --query-cover ${params.dmnd_min_coverage} \
        --id ${params.dmnd_min_identity} \
        --top ${params.dmnd_top_pct} \
        --block-size ${task.memory.toMega() / (1024 * 6 )} \
        --query-gencode ${params.gencode} \
        --compress 1 \
        --unal 0
    done

    cat *aln.gz > ${sample_name}.aln.gz
    """

}


// Filter the alignments with the FAMLI algorithm
process Famli {
    tag "Deduplicate multi-mapping reads"
    container "${container__FAMLI}"
    label 'mem_veryhigh'
    publishDir "${params.output_folder}/abund/details/", mode: "copy"
    errorStrategy 'finish'
    
    input:
    tuple val(sample_name), file(input_aln)
    
    output:
    path "${sample_name}.json.gz"

    """
    set -e
    famli \
      filter \
      --input ${input_aln} \
      --output ${sample_name}.json \
      --threads ${task.cpus} \
      --batchsize ${params.famli_batchsize} \
      --sd-mean-cutoff ${params.sd_mean_cutoff}
    gzip ${sample_name}.json
    """

}

process AssembleAnndata {
    tag "Make gene ~ sample abundance matrix in anndata format"
    container "${container__anndata}"
    label "mem_veryhigh"
    errorStrategy 'finish'
    publishDir "${params.output_folder}/abund/gene/", mode: "copy"
    
    input:
        path sample_jsons
        path peptide_fasta

    output:
        path "specimen_peptide_quantity.h5ad"

"""
#!/usr/bin/env python3

import logging
import numpy as np
import os
import pandas as pd
import json
import gzip
import fastalite
import anndata as ad
from scipy import sparse

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [AssembleAbundances] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# -- Actual work starts here!

jsons = "${sample_jsons.join(";;")}".split(';;')

# First loop through to get specimens and alleles
logging.info(f"Loading in specimens and peptides from {len(jsons):,d} specimens")
alleles = set()
specimens = set()
for fn in jsons:
    specimen = os.path.basename(fn).replace('.json.gz', '')
    specimens.add(specimen)
    sp_json = json.load(
        gzip.open(fn, 'rt')
    )
    alleles.update({
        v['id']
        for v in sp_json
    })
logging.info(f'Specimens {len(specimens):,d}. Peptides {len(alleles):,d}.')
# Convert to list 

specimens = sorted(specimens)
alleles = sorted(alleles)

# and index
specimen_i = {
    v: i for (i, v) in enumerate(specimens)
}
allele_i = {
     v: i for (i, v) in enumerate(alleles)
}

logging.info("Now second loop to extract data")
nreads = []
coverage = []
length = []
depth = []
std = []
i = []
j = []
# Second loop, this to build the coordinates data lists:
# [data, (i, j)]
for fn in jsons:
    specimen = os.path.basename(fn).replace('.json.gz', '')
    sp_i = specimen_i.get(specimen)
    
    sp_json = json.load(
        gzip.open(fn, 'rt')
    )
    i += [sp_i] * len(sp_json)
    j += [
        allele_i.get(v['id'])
        for v in sp_json
    ]
    nreads += [
        v['nreads']
        for v in sp_json
    ]
    coverage += [
        v['coverage']
        for v in sp_json
    ]
    length += [
        v['length']
        for v in sp_json
    ]
    depth += [
        v['depth']
        for v in sp_json
    ]
    std += [
        v['std']
        for v in sp_json
    ]
logging.info("Now generating the anndata object with sparse matrix")
ad_wgs = ad.AnnData(
    sparse.coo_matrix(
        (
            nreads,
            (
                i,
                j
            )
        ),
        shape=(len(specimens), len(alleles)),
        dtype=int
    ).tocsr(),
    obs=pd.DataFrame(index=specimens),
    var=pd.DataFrame(index=alleles),
    dtype=int
)
ad_wgs.layers['coverage'] = sparse.coo_matrix(
    (
        coverage,
        (
            i,
            j
        )
    ),
    shape=(len(specimens), len(alleles)),
    dtype=int
).tocsr()

ad_wgs.layers['length'] = sparse.coo_matrix(
    (
        length,
        (
            i,
            j
        )
    ),
    shape=(len(specimens), len(alleles)),
    dtype=int
).tocsr()

ad_wgs.layers['depth'] = sparse.coo_matrix(
    (
        depth,
        (
            i,
            j
        )
    ),
    shape=(len(specimens), len(alleles)),
    dtype=int
).tocsr()

ad_wgs.layers['std'] = sparse.coo_matrix(
    (
        std,
        (
            i,
            j
        )
    ),
    shape=(len(specimens), len(alleles)),
    dtype=int
).tocsr()

logging.info("Adding in some per-peptide aggregate data")
ad_wgs.layers['detected'] = ad_wgs.X > 0
ad_wgs.var['n_detected'] = np.ravel(ad_wgs.layers['detected'].sum(axis=0))
ad_wgs.var['prevalence'] = ad_wgs.var['n_detected'] / len(specimens)
ad_wgs.var['mean_coverage'] = np.ravel(ad_wgs.layers['coverage'].mean(axis=0))
ad_wgs.var['mean_depth'] = np.ravel(ad_wgs.layers['depth'].mean(axis=0))

logging.info("Injecting peptide sequences into anndata")
peptide_seq = {
    sr.id: sr.seq for sr in fastalite.fastalite(
        fastalite.Opener()('${peptide_fasta}')
    )
}

ad_wgs.var['seq'] = ad_wgs.var.index.map(peptide_seq.get)

logging.info("Writing to H5AD format")
ad_wgs.write_h5ad('specimen_peptide_quantity.h5ad')

logging.info("Completed")

"""

}


//
// Steps to run alignment independently.
//
// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.help = false
params.output = './results/'
params.manifest = null
params.output_prefix = 'geneshot'
params.gene_fasta = null



// Make sure that --output ends with trailing "/" characters
if (!params.output.endsWith("/")){
    params.output_folder = params.output.concat("/")
} else {
    params.output_folder = params.output
}

// imports
include { Read_manifest } from './general'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run Golob-Minot/geneshot/modules/alignment.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)
      --gene_fasta          Compressed FASTA with pre-generated catalog of microbial genes.


    Options:
      --output              Folder to place analysis outputs (default ./results)
      --output_prefix       Text used as a prefix for summary HDF5 output files (default: geneshot)
      -w                    Working directory. Defaults to `./work`

    For Alignment:
      --dmnd_min_identity   Amino acid identity cutoff used to align short reads (default: 90) (DIAMOND)
      --dmnd_min_coverage   Query coverage cutoff used to align short reads (default: 50) (DIAMOND)
      --dmnd_top_pct        Keep top X% of alignments for each short read (default: 1) (DIAMOND)
      --dmnd_min_score      Minimum score for short read alignment (default: 20) (DIAMOND)
      --gencode             Genetic code used for conceptual translation (default: 11) (DIAMOND)
      --sd_mean_cutoff      Ratio of standard deviation / mean depth of sequencing used to filter genes (default: 3.0) (FAMLI)
      --famli_batchsize     Number of alignments to deduplicate in batches (default: 10000000) (FAMLI)

    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This CANNOT be repeated. One row per specimen.
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      These MUST be preprocessed for this workflow.
    """.stripIndent()
}

workflow {
    main:

 
    // Show help message if the user specifies the --help flag at runtime
    if (params.help || params.manifest == null || params.gene_fasta == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }
    // Read and validate manifest
    manifest_file = Channel.from(file(params.manifest))
    manifest = Read_manifest(manifest_file).valid_paired


    Alignment_wf(
        file( params.gene_fasta ),
        manifest.map{ [it.specimen, file(it.R1), file(it.R2)] },
        params.output_prefix
  )

}