// Processes used for alignment of reads against gene databases

workflow alignment_wf {
    get:
        gene_fasta
        reads_ch

    main:

    // Make a DIAMOND indexed database from those gene sequences
    makeDiamondDB(
        gene_fasta
    )

    // Align all specimens against the DIAMOND database
    diamond(
        reads_ch,
        makeDiamondDB.out
    )

    // Filter to the most likely single alignment per query
    famli(
        diamond.out
    )

    // Make a single table with the abundance of every gene across every sample
    assembleAbundances(
        famli.out.toSortedList()
    )

    // Group genes into Co-Abundant Gene Groups (CAGs)
    makeCAGs(
        assembleAbundances.out
    )

    // Calculate the relative abundance of each CAG in these samples
    calcCAGabund(
        assembleAbundances.out,
        makeCAGs.out
    )

    emit:
        cag_csv = makeCAGs.out
        gene_abund_feather = assembleAbundances.out
        cag_abund_feather = calcCAGabund.out
        famli_json_list = famli.out.toSortedList()

}

// Align each sample against the reference database of genes using DIAMOND
process makeDiamondDB {
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

    cat ${R1} ${R2} > query.fastq.gz

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
    container "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
    label "mem_veryhigh"
    errorStrategy 'retry'
    publishDir "${params.output_folder}/abund/", mode: "copy"

    input:
    file sample_jsons

    output:
    file "gene.abund.feather"


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

# Write out to a feather file
logging.info("Writing to disk")
df.reset_index(inplace=True)
df.to_feather("gene.abund.feather")

logging.info("Done")

    """

}

// Make CAGs for each set of samples, at each level of clustering
process makeCAGs {
    container "quay.io/fhcrc-microbiome/find-cags@sha256:30ec0e8b25ef142b68eebcfa84a7a2eeb44ebb25481a60db923fed288abec4a9"
    label "mem_veryhigh"
    errorStrategy 'retry'
    publishDir "${params.output_folder}/ref/", mode: "copy"

    input:
    path gene_feather

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
    '%(asctime)s %(levelname)-8s [makeCAGs] %(message)s'
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

# Set the file path
gene_feather = "${gene_feather}"

# Make sure the file exists
assert os.path.exists(gene_feather), gene_feather

logging.info("Reading in gene abundances from %s" % gene_feather)
df = feather.read_dataframe(
    gene_feather
).set_index("index").applymap(
    np.float16
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