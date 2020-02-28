container__find_cags = "quay.io/fhcrc-microbiome/find-cags:v0.12.1"

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