// containers
container__FAMLI = 'quay.io/fhcrc-microbiome/famli:v1.5'
container__collection = 'quay.io/fhcrc-microbiome/experiment-collection:v0.2'

// Processes used for alignment of reads against gene databases

params.cag_batchsize = 10000

// Default options
params.distance_threshold = 0.5
params.distance_metric = "cosine"
params.linkage_type = "average"
params.famli_batchsize = 10000000

include { makeInitialCAGs } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round1 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round2 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round3 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round4 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round5 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round6 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round7 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round8 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round9 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_round10 } from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include { refineCAGs as refineCAGs_final } from "./make_cags" params(
    distance_threshold: params.distance_threshold,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)

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
}

workflow CAG_contig_oriented_wf {
    //
    //  CAG Stuff starts here
    //    This is for the 'contig oriented'
    //
    take:
        gene_abundances_zarr_tar
        gene_lists

    main: 
    // Group shards of genes into Co-Abundant Gene Groups (CAGs)
    makeInitialCAGs(
        gene_abundances_zarr_tar,
        gene_lists
    )

    // Perform multiple rounds of combining shards to make ever-larger CAGs
    refineCAGs_round1(
        makeInitialCAGs.out[0].toSortedList().flatten().collate(2),
        makeInitialCAGs.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round2(
        refineCAGs_round1.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round1.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round3(
        refineCAGs_round2.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round2.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round4(
        refineCAGs_round3.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round3.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round5(
        refineCAGs_round4.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round4.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round6(
        refineCAGs_round5.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round5.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round7(
        refineCAGs_round6.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round6.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round8(
        refineCAGs_round7.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round7.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round9(
        refineCAGs_round8.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round8.out[1].toSortedList().flatten().collate(2),
    )
    refineCAGs_round10(
        refineCAGs_round9.out[0].toSortedList().flatten().collate(2),
        refineCAGs_round9.out[1].toSortedList().flatten().collate(2),
    )

    // Combine the shards and make a new set of CAGs
    refineCAGs_final(
        refineCAGs_round10.out[0].toSortedList(),
        refineCAGs_round10.out[1].toSortedList(),
    )

    emit:
        cag_csv = refineCAGs_final.out[0]
        cag_abund_feather = refineCAGs_final.out[1]
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


// Make a single feather file with the abundance of every gene across every sample
process AssembleAbundances {
    tag "Make gene ~ sample abundance matrix"
    container "${container__collection}"
    label "mem_veryhigh"
    errorStrategy 'finish'

    input:
    file sample_jsons
    val cag_batchsize
    val output_prefix

    output:
    path "gene_list.*.csv.gz", emit: gene_lists
    path "specimen_gene_count.csv.gz", emit: specimen_gene_count_csv
    path "${output_prefix}.details.hdf5", emit: detailed_hdf
    path "gene_length.csv.gz", emit: gene_length_csv
    path "specimen_reads_aligned.csv.gz", emit: specimen_reads_aligned_csv
    path "gene_abundance.zarr.tar", emit: abundance_zarr_tar


    """
#!/usr/bin/env python3

import logging
import numpy as np
import os
import pandas as pd
import tarfile
import gzip
import json
import pickle
import zarr
pickle.HIGHEST_PROTOCOL = 4

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

def read_json(fp):
    return pd.DataFrame(
        json.load(
            gzip.open(fp, "rt")
        )
    )

# Parse the file list
sample_jsons = "${sample_jsons}".split(" ")
logging.info("Getting ready to read in data for %d samples" % len(sample_jsons))

# Keep track of the complete set of sample names observed in these samples
all_sample_names = set([])

# Also keep track of the complete set of gene names observed in these samples
all_gene_names = set([])

# All abundance tables will be written out to HDF5
store = pd.HDFStore("${output_prefix}.details.hdf5", "w")

# Keep track of the length of each gene
gene_length_dict = dict()

# Keep track of the number of reads aligned per sample
specimen_reads_aligned = dict()

# Keep track of the number of genes detected per sample
specimen_genes_detected = dict()

# Iterate over the list of files
for fp in sample_jsons:
    # Get the sample name from the file name
    # This is only possible because we control the file naming scheme in the `famli` process
    assert fp.endswith(".json.gz")
    sample_name = fp[:-len(".json.gz")]

    logging.info("Reading in %s from %s" % (sample_name, fp))
    df = read_json(fp)

    logging.info("Saving to HDF5")
    df.to_hdf(store, "/abund/gene/long/%s" % sample_name)
    
    # Add the sample name to the total list
    all_sample_names.add(sample_name)

    # Add the gene names to the total list
    for gene_name in df["id"].values:
        all_gene_names.add(gene_name)

    # Add the gene lengths
    for _, r in df.iterrows():
        gene_length_dict[r["id"]] = r["length"]

    # Add in the number of reads aligned
    specimen_reads_aligned[sample_name] = df["nreads"].sum()

    # Add in the number of genes detected
    specimen_genes_detected[sample_name] = df.shape[0]

store.close()

# Serialize sample and gene names
sample_names = list(all_sample_names)
all_gene_names = list(all_gene_names)

# Write out the number of genes detected per sample
pd.DataFrame(
    dict(
        n_genes_aligned = specimen_genes_detected
    )
).reset_index(
).rename(
    columns = dict(
        index = "specimen"
    )
).to_csv(
    "specimen_gene_count.csv.gz",
    index=None,
    compression = "gzip"
)

# Write out the number of reads aligned per sample
pd.DataFrame(
    dict(
        n_reads_aligned = specimen_reads_aligned
    )
).reset_index(
).rename(
    columns = dict(
        index = "specimen"
    )
).to_csv(
    "specimen_reads_aligned.csv.gz",
    index=None,
    compression = "gzip"
)

# Write out the CSV table with the length of each gene
pd.DataFrame([
    dict(
        gene = gene, 
        length = length
    )
    for gene, length in gene_length_dict.items()
]).to_csv(
    "gene_length.csv.gz",
    index = None,
    compression = "gzip"
)

# Write out the gene names in batches of ${cag_batchsize}
for ix, gene_list in enumerate([
    all_gene_names[ix: (ix + ${cag_batchsize})]
    for ix in range(0, len(all_gene_names), ${cag_batchsize})
]):
    print("Writing out %d genes in batch %d" % (len(gene_list), ix))
    with gzip.open("gene_list.%d.csv.gz" % ix, "wt") as handle:
        handle.write("\\n".join(gene_list))

# Write out the sequencing depth / rel abund of each gene in each specimen in zarr format
z_ra = zarr.open(
    "gene_abundance.zarr",
    mode="w",
    shape=(len(all_gene_names), len(sample_names)), 
    chunks=True,
    dtype='f4'
)

z_nr = zarr.open(
    "gene_nreads.zarr",
    mode="w",
    shape=(len(all_gene_names), len(sample_names)), 
    chunks=True,
    dtype='ulonglong'
)

# Iterate over the list of files
for fp in sample_jsons:
    # Get the sample name from the file name
    assert fp.endswith(".json.gz")
    sample_name = fp[:-len(".json.gz")]

    logging.info("Reading in %s from %s" % (sample_name, fp))
    df = read_json(fp)

    # Calculate the depth and rel abund of sequencing per gene
    gene_depth = df.set_index(
        "id"
    ).reindex(
        index=all_gene_names
    )[
        "depth"
    ].fillna(
        0
    )
    gene_ra = gene_depth / df[
        "depth"
    ].sum()

    logging.info("Saving %s to gene_abundance.zarr" % sample_name)
    # Save the sequencing depth to zarr
    z_nr[
        :,
        sample_names.index(sample_name)
    ] = gene_depth.values

    z_ra[
        :,
        sample_names.index(sample_name)
    ] = gene_ra.values

# Write out the sample names and gene names
logging.info("Writing out sample_names.json.gz")
with gzip.open("sample_names.json.gz", "wt") as fo:
    json.dump(sample_names, fo)

logging.info("Writing out gene_names.json.gz")
with gzip.open("gene_names.json.gz", "wt") as fo:
    json.dump(all_gene_names, fo)

logging.info("Creating gene_abundance.zarr.tar")
with tarfile.open("gene_abundance.zarr.tar", "w") as tar:
    for name in [
        "gene_names.json.gz",
        "sample_names.json.gz",
        "gene_abundance.zarr",
        "gene_nreads.zarr",
    ]:
        logging.info("Adding %s to gene_abundance.zarr.tar" % name)
        tar.add(name)

logging.info("Done")

    """

}

process OutputSpecimenGeneAbundance {
    tag "Output specimen-gene-abundances"
    container "${container__collection}"
    label "mem_veryhigh"
    errorStrategy 'finish'
    publishDir "${params.output_folder}/abund/gene/", mode: "copy"

    input:
        path gene_abundances_zarr_tar
    
    output:
        path gene_abundances_zarr_tar
        path 'specimen_gene_abundance_long.csv.gz'
    
"""
#!/usr/bin/env python3
import tarfile
import zarr
import gzip
import csv
import os
import json

with tarfile.open("${gene_abundances_zarr_tar}") as tar:
    tar.extractall()

# Make sure that the expected contents are present
for fp in [
    "gene_names.json.gz",
    "sample_names.json.gz",
    "gene_abundance.zarr",
    "gene_nreads.zarr",
]:
    assert os.path.exists(fp)

# Read in the gene names and sample names as indexed in the zarr
with gzip.open("gene_names.json.gz", "rt") as handle:
    gene_names = json.load(handle)
with gzip.open("sample_names.json.gz", "rt") as handle:
    sample_names = json.load(handle)

z_ra = zarr.open("gene_abundance.zarr", mode="r")
z_rc = zarr.open("gene_nreads.zarr", mode="r")

with gzip.open('specimen_gene_abundance_long.csv.gz', 'wt') as out_h:
    out_w = csv.writer(out_h)
    out_w.writerow(
        (
            'specimen',
            'gene',
            'rel_abd',
            'nreads'
        )
    )
    for gene_ix, gene_name in enumerate(gene_names):
        gene_ra = z_ra[gene_ix, :]
        gene_rc = z_rc[gene_ix, :]
        out_w.writerows((
            (
                sample_names[sam_i],
                gene_name,
                sg_ra,
                sg_rc,
            )
            for sam_i, (sg_ra, sg_rc) in enumerate(zip(gene_ra, gene_rc)) 
            if sg_rc > 0
        ))
"""

}

// Summarize the abundance of every CAG across each sample
process calcCAGabund {
    tag "Make CAG ~ sample abundance matrix"
    container "${container__collection}"
    label "mem_veryhigh"
    errorStrategy 'finish'
    publishDir "${params.output_folder}/abund/", mode: "copy"

    input:
    path cag_csv_gz

    output:
    file "CAG.abund.feather"

    """
#!/usr/bin/env python3

from collections import defaultdict
import gzip
import json
import numpy as np
import os
import pandas as pd
import logging
import tarfile
import zarr

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [calculateCAGabundance] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Read in the dictionary linking genes and CAGs
cags = {
    cag_id: cag_df["gene"].tolist()
    for cag_id, cag_df in pd.read_csv(
        "${cag_csv_gz}",
        compression="gzip"
    ).groupby(
        "CAG"
    )
}

# Extract everything from the gene_abundance.zarr.tar
logging.info("Extracting ${gene_abundances_zarr_tar}")
with tarfile.open("${gene_abundances_zarr_tar}") as tar:
    tar.extractall()

# Make sure that the expected contents are present
for fp in [
    "gene_names.json.gz",
    "sample_names.json.gz",
    "gene_abundance.zarr",
]:
    logging.info("Checking that %s is present" % fp)
    assert os.path.exists(fp)

# Read in the gene names and sample names as indexed in the zarr
logging.info("Reading in gene_names.json.gz")
with gzip.open("gene_names.json.gz", "rt") as handle:
    gene_names = json.load(handle)

logging.info("Reading in sample_names.json.gz")
with gzip.open("sample_names.json.gz", "rt") as handle:
    sample_names = json.load(handle)

# Open the zarr store
logging.info("Reading in gene abundances from gene_abundance.zarr")
z = zarr.open("gene_abundance.zarr", mode="r")

# Set up an array for CAG abundances
df = pd.DataFrame(
    data = np.zeros(
        (len(cags), len(sample_names)),
        dtype = np.float32,
    ),
    dtype = np.float32,
    columns = sample_names,
    index = list(range(len(cags))),
)

# Iterate over each sample
for sample_ix, sample_name in enumerate(sample_names):

    # Read the depth of sequencing for each gene
    logging.info("Reading gene abundances for %s" % sample_name)
    sample_gene_depth = pd.Series(z[:, sample_ix], index=gene_names)

    # Sum up the depth by CAG
    df[sample_name] = pd.Series({
        cag_ix: sample_gene_depth.reindex(index=cag_gene_list).sum()
        for cag_ix, cag_gene_list in cags.items()
    })


# Now write out to a feather file
logging.info("Building a single DataFrame of CAG abundances")
df.reset_index(
).to_feather(
    "CAG.abund.feather"
)

logging.info("Done")
    """
}
