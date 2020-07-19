// Processes used for alignment of reads against gene databases

params.cag_batchsize = 100000

// Default options
params.distance_threshold = 0.5
params.distance_metric = "cosine"
params.linkage_type = "average"
params.famli_batchsize = 10000000

include makeInitialCAGs from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include refineCAGs as refineCAGs_round1 from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include refineCAGs as refineCAGs_round2 from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include refineCAGs as refineCAGs_round3 from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include refineCAGs as refineCAGs_round4 from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include refineCAGs as refineCAGs_round5 from "./make_cags" params(
    distance_threshold: params.distance_threshold / 2,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)
include refineCAGs as refineCAGs_final from "./make_cags" params(
    distance_threshold: params.distance_threshold,
    distance_metric: params.distance_metric,
    linkage_type: params.linkage_type
)

workflow alignment_wf {
    take:
        gene_fasta
        reads_ch
        output_prefix

    main:

    // Make a DIAMOND indexed database from those gene sequences
    diamondDB(
        gene_fasta
    )

    // Align all specimens against the DIAMOND database
    diamond(
        reads_ch,
        diamondDB.out
    )

    // Filter to the most likely single alignment per query
    famli(
        diamond.out
    )

    // Make a single table with the abundance of every gene across every sample
    assembleAbundances(
        famli.out.toSortedList(),
        params.cag_batchsize,
        output_prefix
    )

    // Group shards of genes into Co-Abundant Gene Groups (CAGs)
    makeInitialCAGs(
        assembleAbundances.out[2],
        assembleAbundances.out[0].flatten()
    )

    // Perform multiple rounds of combining shards to make ever-larger CAGs
    refineCAGs_round1(
        assembleAbundances.out[2],
        makeInitialCAGs.out.toSortedList().flatten().collate(2)
    )
    refineCAGs_round2(
        assembleAbundances.out[2],
        refineCAGs_round1.out.toSortedList().flatten().collate(2)
    )
    refineCAGs_round3(
        assembleAbundances.out[2],
        refineCAGs_round2.out.toSortedList().flatten().collate(2)
    )
    refineCAGs_round4(
        assembleAbundances.out[2],
        refineCAGs_round3.out.toSortedList().flatten().collate(2)
    )
    refineCAGs_round5(
        assembleAbundances.out[2],
        refineCAGs_round4.out.toSortedList().flatten().collate(2)
    )

    // Combine the shards and make a new set of CAGs
    refineCAGs_final(
        assembleAbundances.out[2],
        refineCAGs_round5.out.collect()
    )

    // Calculate the relative abundance of each CAG in these samples
    calcCAGabund(
        assembleAbundances.out[2],
        refineCAGs_final.out
    )

    emit:
        cag_csv = refineCAGs_final.out
        cag_abund_feather = calcCAGabund.out
        famli_json_list = famli.out.toSortedList()
        specimen_gene_count_csv = assembleAbundances.out[1]
        specimen_reads_aligned_csv = assembleAbundances.out[4]
        detailed_hdf = assembleAbundances.out[2]
        gene_length_csv = assembleAbundances.out[3]
}

// Align each sample against the reference database of genes using DIAMOND
process diamondDB {
    tag "Make a DIAMOND database"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'
    errorStrategy 'retry'
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
process diamond {
    tag "Align to the gene catalog"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    tuple val(sample_name), file(R1), file(R2)
    file refdb
    
    output:
    tuple sample_name, file("${sample_name}.aln.gz")

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
        --block-size ${task.memory.toMega() / (1024 * 6 * task.attempt)} \
        --query-gencode ${params.gencode} \
        --compress 1 \
        --unal 0
    done

    cat *aln.gz > ${sample_name}.aln.gz
    """

}


// Filter the alignments with the FAMLI algorithm
process famli {
    tag "Deduplicate multi-mapping reads"
    container "quay.io/fhcrc-microbiome/famli:v1.5"
    label 'mem_veryhigh'
    publishDir "${params.output_folder}/abund/details/", mode: "copy"
    errorStrategy 'retry'
    
    input:
    tuple sample_name, file(input_aln)
    
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
process assembleAbundances {
    tag "Make gene ~ sample abundance matrix"
    container "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
    label "mem_veryhigh"
    errorStrategy 'retry'

    input:
    file sample_jsons
    val cag_batchsize
    val output_prefix

    output:
    file "gene_list.*.csv.gz"
    file "specimen_gene_count.csv.gz"
    file "${output_prefix}.details.hdf5"
    path "gene_length.csv.gz"
    path "specimen_reads_aligned.csv.gz"


    """
#!/usr/bin/env python3

import logging
import numpy as np
import os
import pandas as pd
import gzip
import json
import pickle
pickle.HIGHEST_PROTOCOL = 4

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

# Write out the gene names in batches of ${params.cag_batchsize}
for ix, gene_list in enumerate([
    all_gene_names[ix: (ix + ${cag_batchsize})]
    for ix in range(0, len(all_gene_names), ${cag_batchsize})
]):
    print("Writing out %d genes in batch %d" % (len(gene_list), ix))
    with gzip.open("gene_list.%d.csv.gz" % ix, "wt") as handle:
        handle.write("\\n".join(gene_list))

logging.info("Done")

    """

}


// Summarize the abundance of every CAG across each sample
process calcCAGabund {
    tag "Make CAG ~ sample abundance matrix"
    container "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
    label "mem_veryhigh"
    errorStrategy 'retry'
    publishDir "${params.output_folder}/abund/", mode: "copy"

    input:
    path details_hdf
    path cag_csv_gz

    output:
    file "CAG.abund.feather"

    """
#!/usr/bin/env python3

from collections import defaultdict
import os
import pandas as pd
import logging

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
cags_dict = pd.read_csv(
    "${cag_csv_gz}",
    compression="gzip"
).set_index(
    "gene"
)[
    "CAG"
]

# Set the file path to the HDF5 with gene abundances
details_hdf = "${details_hdf}"

# Make sure the files exist
assert os.path.exists(details_hdf), details_hdf

# Read in the table of gene abundances
abund_df = defaultdict(lambda: defaultdict(float))

# Open the HDF5 store
with pd.HDFStore(details_hdf, "r") as store:

    # Iterate over keys in the store
    for k in store:

        # Sample abundances have a certain prefix
        if k.startswith("/abund/gene/long/"):

            # Parse the sample name
            sample_name = k.replace("/abund/gene/long/", "")
            logging.info("Reading gene abundances for {}".format(sample_name))

            # Read the depth of sequencing for each gene
            sample_depth = pd.read_hdf(
                store,
                k
            ).set_index(
                "id"
            )[
                "depth"
            ]

            # Normalize to sum to 1
            sample_depth = sample_depth / sample_depth.sum()

            # Add the abundance of each gene to its linked CAG
            for gene_name, gene_prop in sample_depth.items():
                abund_df[
                    sample_name
                ][
                    cags_dict[
                        gene_name
                    ]
                ] += gene_prop

# Now write out to a feather file
logging.info("Building a single DataFrame of CAG abundances")
pd.DataFrame(
    abund_df
).fillna(
    0
).reset_index(
).to_feather(
    "CAG.abund.feather"
)

logging.info("Done")
    """
}