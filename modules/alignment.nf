// Processes used for alignment of reads against gene databases

params.cag_batchsize = 1000000

workflow alignment_wf {
    take:
        gene_fasta
        reads_ch

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
        params.cag_batchsize
    )

    // Group shards of genes into Co-Abundant Gene Groups (CAGs)
    makeInitialCAGs(
        assembleAbundances.out[0],
        assembleAbundances.out[1].flatten()
    )

    // Combine the shards and make a new set of CAGs
    makeFinalCAGs(
        assembleAbundances.out[0],
        makeInitialCAGs.out.collect()
    )

    // Calculate the relative abundance of each CAG in these samples
    calcCAGabund(
        assembleAbundances.out[0],
        makeFinalCAGs.out
    )

    emit:
        cag_csv = makeFinalCAGs.out
        gene_abund_feather = assembleAbundances.out[0]
        cag_abund_feather = calcCAGabund.out
        famli_json_list = famli.out.toSortedList()

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

    cat ${R1} ${R2} | \
    gunzip -c | \
    awk "{if(NR % 4 == 1){print \\"@\\" NR }else{if(NR % 4 == 3){print \\"+\\"}else{print}}}" | \
    gzip -c \
     > query.fastq.gz

    diamond \
      blastx \
      --query query.fastq.gz \
      --out ${sample_name}.aln.gz \
      --threads ${task.cpus} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --min-score ${params.dmnd_min_score} \
      --query-cover ${params.dmnd_min_coverage} \
      --id ${params.dmnd_min_identity} \
      --top ${params.dmnd_top_pct} \
      --block-size ${task.memory.toMega() / (1024 * 6)} \
      --query-gencode ${params.gencode} \
      --compress 1 \
      --unal 0
    """

}


// Filter the alignments with the FAMLI algorithm
process famli {
    tag "Deduplicate multi-mapping reads"
    container "quay.io/fhcrc-microbiome/famli@sha256:241a7db60cb735abd59f4829e8ddda0451622b6eb2321f176fd9d76297d8c9e7"
    label 'mem_veryhigh'
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
      --batchsize 5000000
    gzip ${sample_name}.json
    """

}


// Make a single feather file with the abundance of every gene across every sample
process assembleAbundances {
    tag "Make gene ~ sample abundance matrix"
    container "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
    label "mem_veryhigh"
    errorStrategy 'retry'
    publishDir "${params.output_folder}/abund/", mode: "copy"

    input:
    file sample_jsons
    val cag_batchsize

    output:
    file "gene.abund.feather"
    file "gene_list.*.csv.gz"


    """
#!/usr/bin/env python3

import logging
import numpy as np
import os
import pandas as pd
import gzip
import json

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
    return {
        r["id"]: r["depth"]
        for r in json.load(gzip.open(fp, "rt"))
    }

# Parse the file list
sample_jsons = "${sample_jsons}".split(" ")
logging.info("Getting ready to read in data for %d sample" % len(sample_jsons))

# All of the abundances will go into a single dict
all_abund = dict()
# Also keep track of the complete set of gene names observed in these samples
all_gene_names = set([])

# Iterate over the list of files
for fp in sample_jsons:
    # Get the sample name from the file name
    # This is only possible because we control the file naming scheme in the `famli` process
    assert fp.endswith(".json.gz")
    sample_name = fp[:-len(".json.gz")]

    logging.info("Reading in %s from %s" % (sample_name, fp))
    all_abund[sample_name] = read_json(fp)

    # Add the gene names to the total list
    for gene_name in all_abund[sample_name]:
        all_gene_names.add(gene_name)

# Serialize sample and gene names
sample_names = list(all_abund.keys())
all_gene_names = list(all_gene_names)

# Now make a DataFrame for all of this data
df = pd.DataFrame(
    np.zeros((len(all_gene_names), len(sample_names)), dtype=np.float32),
    columns=sample_names,
    index=all_gene_names
)

# Add all of the data to the DataFrame
for sample_name, sample_abund in all_abund.items():
    # Format as a Pandas Series
    sample_abund = pd.Series(
        sample_abund
    )

    # Divide by the sum to get the proportional abundance
    sample_abund = sample_abund / sample_abund.sum()

    # Add the values to the table
    df.values[
        :, sample_names.index(sample_name)
    ] = sample_abund.reindex(
        index=all_gene_names
    ).fillna(
        0
    ).apply(
        np.float32
    ).values

# Write out the abundances to a feather file
logging.info("Writing to disk")
df.reset_index(inplace=True)
df.to_feather("gene.abund.feather")

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

// Make CAGs for each set of samples, with the subset of genes for this shard
process makeInitialCAGs {
    tag "Group gene subsets by co-abundance"
    container "quay.io/fhcrc-microbiome/find-cags@sha256:30ec0e8b25ef142b68eebcfa84a7a2eeb44ebb25481a60db923fed288abec4a9"
    label "mem_veryhigh"
    errorStrategy 'retry'

    input:
    path gene_feather
    path gene_list_csv

    output:
    file "CAGs.csv.gz"

    """
#!/usr/bin/env python3

import feather
import gzip
import json
import logging
from multiprocessing import Pool
import nmslib
import numpy as np
import os
import pandas as pd
from ann_linkage_clustering.lib import make_cags_with_ann
from ann_linkage_clustering.lib import iteratively_refine_cags
from ann_linkage_clustering.lib import make_nmslib_index

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [makeInitialCAGs] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Set up the multiprocessing pool
threads = int("${task.cpus}")
pool = Pool(threads)

# Set the file path to the feather with gene abundances
gene_feather = "${gene_feather}"

# Set the file path to the genes for this subset
gene_list_csv = "${gene_list_csv}"

# Make sure the files exist
assert os.path.exists(gene_feather), gene_feather
assert os.path.exists(gene_list_csv), gene_list_csv

logging.info("Reading in gene abundances from %s" % gene_feather)
df = feather.read_dataframe(
    gene_feather
).set_index("index").applymap(
    np.float16
)

logging.info("Reading in the list of genes for this shard from %s" % (gene_list_csv))
gene_list = [
    line.rstrip("\\n")
    for line in gzip.open(gene_list_csv, "rt")
]
logging.info("This shard contains %d genes" % (len(gene_list)))

# Subset the gene abundances
logging.info("Making sure that the gene names match between the feather and CSV")
n_mismatch = len(set(gene_list) - set(df.index.values))
msg = "Gene names do not match between the feather and CSV (n= %d / %d)" % (n_mismatch, len(gene_list))
assert n_mismatch == 0, msg
df = df.reindex(index=gene_list)
logging.info("Subset down to %d genes for this shard" % (len(gene_list)))

max_dist = float("${params.distance_threshold}")
logging.info("Maximum cosine distance: %s" % max_dist)

# Make the nmslib index
logging.info("Making the HNSW index")
index = nmslib.init(method='hnsw', space='cosinesimil')
logging.info("Adding %d genes to the nmslib index" % (df.shape[0]))
index.addDataPointBatch(df.values)
logging.info("Making the index")
index.createIndex({'post': 2, "M": 100}, print_progress=True)


# Make the CAGs
logging.info("Making first-round CAGs")
cags = make_cags_with_ann(
    index,
    max_dist,
    df.copy(),
    pool,
    threads=threads,
    distance_metric="${params.distance_metric}",
    linkage_type="${params.linkage_type}"
)

logging.info("Closing the process pool")
pool.close()

logging.info("Clearing the previous index from memory")
del index

logging.info("Refining CAGS")
iteratively_refine_cags(
    cags,
    df.copy(),
    max_dist,
    threads=threads,
    distance_metric="${params.distance_metric}",
    linkage_type="${params.linkage_type}"
)

logging.info("Formatting CAGs as a DataFrame")
cags_df = pd.DataFrame(
    [
        [ix, gene_id]
        for ix, list_of_genes in enumerate(
            sorted(
                list(
                    cags.values()
                ), 
                key=len, 
                reverse=True
            )
        )
        for gene_id in list_of_genes
    ],
    columns=["CAG", "gene"]
)

# Make sure that we have every gene assigned to a CAG
assert cags_df.shape[0] == df.shape[0], (cags_df.shape[0], df.shape[0])

logging.info("Largest CAGs:")
print(cags_df["CAG"].value_counts().head())

fp_out = "CAGs.csv.gz"

logging.info("Writing out CAGs to %s" % fp_out)
cags_df.to_csv(fp_out, compression="gzip", index=None)

logging.info("Done")
    """
}

// Make CAGs for each set of samples, combining CAGs made for each individual shard
process makeFinalCAGs {
    tag "Group all genes by co-abundance"
    container "quay.io/fhcrc-microbiome/find-cags@sha256:30ec0e8b25ef142b68eebcfa84a7a2eeb44ebb25481a60db923fed288abec4a9"
    label "mem_veryhigh"
    errorStrategy 'retry'
    publishDir "${params.output_folder}/ref/", mode: "copy"

    input:
    path gene_feather
    path "shard.CAG.*.csv.gz"

    output:
    file "CAGs.csv.gz"

    """
#!/usr/bin/env python3

import feather
import gzip
import json
import logging
from multiprocessing import Pool
import nmslib
import numpy as np
import os
import pandas as pd
from ann_linkage_clustering.lib import make_cags_with_ann
from ann_linkage_clustering.lib import iteratively_refine_cags
from ann_linkage_clustering.lib import make_nmslib_index

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [makeFinalCAGs] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Set up the multiprocessing pool
threads = int("${task.cpus}")
pool = Pool(threads)

# Set the file path to the feather with gene abundances
gene_feather = "${gene_feather}"

# Make sure the files exist
assert os.path.exists(gene_feather), gene_feather

# Set the file path to the CAGs made for each subset
cag_csv_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("shard.CAG.")
]
# Make sure all of the files have the complete ending
# (incompletely staged files will have a random suffix appended)
for fp in cag_csv_list:
    assert fp.endswith(".csv.gz"), "Incomplete input file found: %s" % (fp)

assert len(cag_csv_list) > 0, "Didn't find CAGs from any previous shard"

logging.info("Reading in gene abundances from %s" % gene_feather)
df = feather.read_dataframe(
    gene_feather
).set_index("index").applymap(
    np.float16
)

logging.info("Reading in CAGs from previous shard")
cags = dict()
ix = 0
n = 0
for fp in cag_csv_list:
    shard_cags = pd.read_csv(fp, compression="gzip", sep=",")
    for _, gene_list in shard_cags.groupby("CAG"):
        cags[ix] = gene_list['gene'].tolist()
        ix += 1
        n += gene_list.shape[0]
logging.info("Read in %d CAGs from %d shards covering %d genes" % (len(cags), len(cag_csv_list), n))

max_dist = float("${params.distance_threshold}")
logging.info("Maximum cosine distance: %s" % max_dist)

logging.info("Refining CAGS")
iteratively_refine_cags(
    cags,
    df.copy(),
    max_dist,
    threads=threads,
    distance_metric="${params.distance_metric}",
    linkage_type="${params.linkage_type}"
)

logging.info("Formatting CAGs as a DataFrame")
cags_df = pd.DataFrame(
    [
        [ix, gene_id]
        for ix, list_of_genes in enumerate(
            sorted(
                list(
                    cags.values()
                ), 
                key=len, 
                reverse=True
            )
        )
        for gene_id in list_of_genes
    ],
    columns=["CAG", "gene"]
)

# Make sure that we have every gene assigned to a CAG
assert cags_df.shape[0] == df.shape[0], (cags_df.shape[0], df.shape[0])

logging.info("Largest CAGs:")
print(cags_df["CAG"].value_counts().head())

fp_out = "CAGs.csv.gz"

logging.info("Writing out CAGs to %s" % fp_out)
cags_df.to_csv(fp_out, compression="gzip", index=None)

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
    path gene_feather
    path cag_csv_gz

    output:
    file "CAG.abund.feather"

    """
#!/usr/bin/env python3

import feather
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

# Read in the table of CAGs
cags_df = pd.read_csv(
    "${cag_csv_gz}",
    compression="gzip"
)

# Read in the table of gene abundances
abund_df = pd.read_feather(
    "${gene_feather}"
).set_index(
    "index"
)

# Annotate each gene with the CAG it was assigned to
abund_df["CAG"] = cags_df.set_index("gene")["CAG"]

# Make sure that every gene was assigned to a CAG
assert abund_df["CAG"].isnull().sum() == 0

# Now sum up the gene relative abundance by CAGs
# and write out to a feather file
abund_df.groupby(
    "CAG"
).sum(
).reset_index(
).to_feather(
    "CAG.abund.feather"
)

logging.info("Done")
    """
}