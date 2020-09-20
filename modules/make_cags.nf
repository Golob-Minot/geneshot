container__find_cags = "quay.io/fhcrc-microbiome/find-cags:v0.13.0"

// Default options
params.distance_threshold = 0.5
params.distance_metric = "cosine"
params.linkage_type = "average"

// Make CAGs for each set of samples, with the subset of genes for this shard
process makeInitialCAGs {
    tag "Group gene subsets by co-abundance"
    container "${container__find_cags}"
    label "mem_medium"
    errorStrategy 'retry'

    input:
    path gene_abundances_zarr_tar
    path gene_list_csv

    output:
    file "CAGs.csv.gz"
    file "CAGs.abund.feather"

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
import shutil
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
    data = np.zeros(
        (len(gene_list), len(sample_names)),
        dtype = np.float32,
    ),
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

logging.info("Sorting CAGs by size")
cags = {
    ix: list_of_genes
    for ix, list_of_genes in enumerate(
        sorted(
            list(
                cags.values()
            ), 
            key=len, 
            reverse=True
        )
    )
}

logging.info("Computing relative abundance of CAGs")
# Rows are CAGs, columns are specimens
cags_abund_df = pd.DataFrame({
    cag_id: df.reindex(
        index=list_of_genes
    ).sum()
    for cag_id, list_of_genes in cags.items()
}).T.sort_index()

fp_out = "CAGs.abund.feather"
logging.info("Saving CAG relative abundances to %s" % fp_out)
cags_abund_df.reset_index().to_feather(
    fp_out
)

logging.info("Formatting CAG membership as a DataFrame")
cags_df = pd.DataFrame(
    [
        [ix, gene_id]
        for ix, list_of_genes in cags.items()
        for gene_id in list_of_genes
    ],
    columns=["CAG", "gene"]
)

logging.info("Largest CAGs:")
for cag_id, cag_size in cags_df["CAG"].value_counts().head().items():
    logging.info("CAG ID: %d, %d genes" % (cag_id, cag_size))

fp_out = "CAGs.csv.gz"

logging.info("Writing out CAG membership to %s" % fp_out)
cags_df.to_csv(fp_out, compression="gzip", index=None)

logging.info("Deleting the temporary zarr")
del z
shutil.rmtree("gene_abundance.zarr")

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
    file "shard.CAG.*.csv.gz"
    file "shard.CAG.*.feather"

    output:
    file "CAGs.csv.gz"
    file "CAGs.abund.feather"

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
import tarfile
import zarr
import shutil
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

# Set the file path to the CAGs made for each subset
cag_csv_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("shard.CAG.") and ".csv.gz" in fp
]
cag_feather_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("shard.CAG.") and ".feather" in fp
]

# Make sure all of the files have the complete ending
# (incompletely staged files will have a random suffix appended)
for fp in cag_csv_list:
    assert fp.endswith(".csv.gz"), "Incomplete input file found: %s" % (fp)

    # Make sure that we have a feather file which matches the CAG membership
    assert fp.replace(".csv.gz", ".feather") in cag_feather_list, "Names of CAG membership and feather files do not match"

assert len(cag_csv_list) > 0, "Didn't find CAGs from any previous shard"
assert len(cag_csv_list) == len(cag_feather_list), "Number of CSV and feather files does not match"

# Keep track of the abundances of the CAGs from the inputs
cag_abund = []

# Keep track of the genes that the previous set of CAGs corresponded to
cag_membership = {}

logging.info("Reading in CAGs from previous shard")
ix = 0
for fp in cag_csv_list:

    # Read in the gene membership from this shard
    shard_cags_membership = pd.read_csv(fp, compression="gzip", sep=",")
    logging.info("Read in %d genes in %d CAGs from %s" % (
        shard_cags_membership.shape[0],
        shard_cags_membership["CAG"].unique().shape[0],
        fp
    ))

    # Read in the abundances of the CAGs in this shard
    feather_fp = fp.replace(".csv.gz", ".feather")
    shard_cags_abundance = pd.read_feather(feather_fp)
    logging.info("Read in abundances for %d CAGs from %s" % (
        shard_cags_abundance.shape[0],
        feather_fp
    ))

    # Make sure that the number of CAGs is the same in both files
    assert shard_cags_abundance.shape[0] == shard_cags_membership["CAG"].unique().shape[0]

    # Transform each of the CAG IDs from the shard into a new CAG ID for the combined set
    cag_id_mapping = {}
    for previous_cag_id in shard_cags_membership["CAG"].unique():
        cag_id_mapping[previous_cag_id] = ix
        ix += 1

    # Record which gene goes with which CAG (with the new IDs)
    for _, r in shard_cags_membership.iterrows():
        cag_membership[r["gene"]] = cag_id_mapping[r["CAG"]]

    # Also change the CAG IDs for the abundance table
    shard_cags_abundance = shard_cags_abundance.replace(
        to_replace={
            "index": cag_id_mapping
        },
    )

    # Add the abundance to the running total
    cag_abund.append(shard_cags_abundance)

# Combine all of the tables
logging.info("Combining CAG abundance across all CAGs")
cag_abund = pd.concat(cag_abund)
cag_membership = pd.Series(cag_membership)

# Calculate the size of all of the input CAGs
input_cag_size = cag_membership.value_counts()

# Make sure that things all add up
logging.info("Number of CAGs with abundances: %d" % cag_abund.shape[0])
logging.info("Number of CAGs with members: %d" % input_cag_size.shape[0])
logging.info("Largest CAG index: %d" % max(input_cag_size.index.values))
assert cag_abund.shape[0] == ix
assert cag_abund.shape[0] == int(max(input_cag_size.index.values))
assert cag_abund.shape[0] == input_cag_size.shape[0]

logging.info(
    "Read in %d CAGs from %d shards covering %d genes" % (
        cag_abund.shape[0], 
        len(cag_csv_list), 
        cag_membership.shape[0]
    )
)


max_dist = float("${params.distance_threshold}")
logging.info("Maximum cosine distance: %s" % max_dist)

# In the `iteratively_refine_cags` step, CAGs will be combined
grouped_cags = {
    cag_ix: [cag_ix]
    for cag_ix in cag_abund.index.values
}

# Parse the number of threads available
threads = int("${task.cpus}")
logging.info("Number of threads available: %d" % threads)

logging.info("Refining CAGS")
iteratively_refine_cags(
    grouped_cags,
    cag_abund.copy(),
    max_dist,
    threads=threads,
    distance_metric="${params.distance_metric}",
    linkage_type="${params.linkage_type}",
    max_iters = 10
)

logging.info("Sorting CAGs by size")
output_cag_size = pd.Series({
    new_cag_id: sum([
        input_cag_size[old_cag_id]
        for old_cag_id in old_cag_id_list
    ])
    for new_cag_id, old_cag_id_list in grouped_cags.items()
}).sort_values(
    ascending=False
)

output_cag_ranking = pd.Series(
    range(output_cag_size.shape[0]), 
    index=output_cag_size.index
)

# Name the CAGs based on the rank order of aggregate size (num. genes)
logging.info("Renaming genes with new CAG groupings")
new_cag_mapping = {
    old_cag_id: output_cag_ranking[new_cag_id]
    for new_cag_id, old_cag_id_list in grouped_cags.items()
    for old_cag_id in old_cag_id_list
}

# Update the CAG membership table and format as a DataFrame
logging.info("Updating the CAG membership table")
cag_membership = pd.DataFrame({
    "gene": cag_membership.index.values,
    "CAG": cag_membership.apply(new_cag_mapping.get)
})

logging.info("Computing the abundance of new CAGs")
cag_abund = pd.DataFrame({
    output_cag_ranking[new_cag_id]: cag_abund.reindex(index=old_cag_id_list).sum()
    for new_cag_id, old_cag_id_list in grouped_cags.items()
}).T

logging.info("Largest CAGs:")
for cag_id, cag_size in cag_membership["CAG"].value_counts().head().items():
    logging.info("CAG %d, %d genes" % (cag_id, cag_size))

fp_out = "CAGs.csv.gz"
logging.info("Writing out CAG membership to %s" % fp_out)
cag_membership.to_csv(fp_out, compression="gzip", index=None)

fp_out = "CAGs.abund.feather"
logging.info("Writing out CAG abundance to %s" % fp_out)
cag_abund.reset_index().to_feather(fp_out)

logging.info("Done")
os._exit(0)
    """
}
