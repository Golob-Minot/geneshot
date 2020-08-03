container__find_cags = "quay.io/fhcrc-microbiome/find-cags:v0.13.0"

// Default options
params.distance_threshold = 0.5
params.distance_metric = "cosine"
params.linkage_type = "average"

// Make CAGs for each set of samples, with the subset of genes for this shard
process makeInitialCAGs {
    tag "Group gene subsets by co-abundance"
    container "${container__find_cags}"
    label "mem_veryhigh"
    errorStrategy 'retry'

    input:
    path gene_abundances_zarr_tar
    path gene_list_csv

    output:
    file "CAGs.csv.gz"

    """
#!/usr/bin/env python3

import gzip
import json
import logging
from multiprocessing import Pool
import nmslib
import numpy as np
import os
import pandas as pd
import tarfile
import zarr
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

# Set the file path to the genes for this subset
gene_list_csv = "${gene_list_csv}"

# Make sure the files exist
assert os.path.exists(details_hdf), details_hdf
assert os.path.exists(gene_list_csv), gene_list_csv

logging.info("Reading in the list of genes for this shard from %s" % (gene_list_csv))
gene_list = [
    line.rstrip("\\n")
    for line in gzip.open(gene_list_csv, "rt")
]
logging.info("This shard contains %d genes" % (len(gene_list)))

# Open the zarr store
logging.info("Reading in gene abundances from gene_abundance.zarr")
z = zarr.open("gene_abundance.zarr", mode="r")

# Set up an array for gene abundances
df = pd.DataFrame(
    shape = (len(gene_list), len(sample_names)),
    dtype = np.float32,
    index = gene_list,
    columns = sample_names
)

# Iterate over each sample
for sample_ix, sample_name in enumerate(sample_names):

    # Read in the abundances for this sample
    logging.info("Reading in gene abundances for %s" % sample_name)
    df[sample_name] = pd.Series(
        z[:, sample_ix],
        index=gene_names
    ).reindex(
        gene_list
    ).apply(
        np.float32
    )

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
    linkage_type="${params.linkage_type}",
    max_iters = 5
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

logging.info("Largest CAGs:")
print(cags_df["CAG"].value_counts().head())

fp_out = "CAGs.csv.gz"

logging.info("Writing out CAGs to %s" % fp_out)
cags_df.to_csv(fp_out, compression="gzip", index=None)

logging.info("Done")
    """
}

// Make CAGs for each set of samples, combining CAGs made for each individual shard
process refineCAGs {
    tag "Group all genes by co-abundance"
    container "${container__find_cags}"
    label "mem_veryhigh"
    errorStrategy 'retry'

    input:
    path gene_abundances_zarr_tar
    path "shard.CAG.*.csv.gz"

    output:
    file "CAGs.csv.gz"

    """
#!/usr/bin/env python3

from collections import defaultdict
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

logging.info("Reading in CAGs from previous shard")
cags = dict()
ix = 0
n = 0
gene_list = set()
for fp in cag_csv_list:
    shard_cags = pd.read_csv(fp, compression="gzip", sep=",")
    for _, shard_df in shard_cags.groupby("CAG"):
        cags[ix] = shard_df['gene'].tolist()
        gene_list.update(set(cags[ix]))
        ix += 1

logging.info(
    "Read in %d CAGs from %d shards covering %d genes" % (
        len(cags), 
        len(cag_csv_list), 
        len(gene_list)
    )
)

# Open the zarr store
logging.info("Reading in gene abundances from gene_abundance.zarr")
z = zarr.open("gene_abundance.zarr", mode="r")

# Set up an array for CAG abundances
df = pd.DataFrame(
    shape = (len(cags), len(sample_names)),
    dtype = np.float32,
    columns = sample_names
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

max_dist = float("${params.distance_threshold}")
logging.info("Maximum cosine distance: %s" % max_dist)

# In the `iteratively_refine_cags` step, CAGs will be combined
grouped_cags = {
    cag_ix: [cag_ix]
    for cag_ix in range(len(cags))
}

logging.info("Refining CAGS")
iteratively_refine_cags(
    grouped_cags,
    df.copy(),
    max_dist,
    threads=threads,
    distance_metric="${params.distance_metric}",
    linkage_type="${params.linkage_type}",
    max_iters = 5
)

logging.info("Expanding with the original set of CAGs")
new_cags = {
    cag_group_ix: [
        gene_name
        for original_cag in cag_group_list
        for gene_name in cags[original_cag]
    ]
    for cag_group_ix, cag_group_list in grouped_cags.items()
}

logging.info("Formatting CAGs as a DataFrame")
cags_df = pd.DataFrame(
    [
        [ix, gene_id]
        for ix, list_of_genes in enumerate(
            sorted(
                list(
                    new_cags.values()
                ), 
                key=len, 
                reverse=True
            )
        )
        for gene_id in list_of_genes
    ],
    columns=["CAG", "gene"]
)

logging.info("Largest CAGs:")
print(cags_df["CAG"].value_counts().head())

fp_out = "CAGs.csv.gz"

logging.info("Writing out CAGs to %s" % fp_out)
cags_df.to_csv(fp_out, compression="gzip", index=None)

logging.info("Done")
    """
}
