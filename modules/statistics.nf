// // Processes to perform statistical analysis

// Container versions
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.2.1_pyarrow"

// Workflow to test out the provided formula and manifest CSV in a dry run
// This is intended to catch the case early on where the provided formula
// is not formatted correctly, or does not match the manifest (sample sheet)
workflow validation_wf {
    take: 
        manifest_csv

    main: 

    // Set up a channel with the strings of the formula(s) provided
    formula_ch = Channel.of(
        params.formula.split(",")
    )

    // Generate mock data
    mockData(
        manifest_csv
    )

    // Split up the mock data to be analyzed in parallel
    splitCorncob(
        mockData.out
    )

    // Run corncob on the mock data
    runCorncob(
        splitCorncob.out.flatten(),
        manifest_csv,
        "CAG",
        formula_ch
    )

    // Join the results
    joinCorncob(
        runCorncob.out.toSortedList(),
        "CAG"
    )

    // Make sure that the results conform to the expected format
    validateFormula(
        joinCorncob.out,
        manifest_csv
    )

    // Return the validated manifest
    emit:
    validateFormula.out

}

// Workflow to run corncob on the actual data
workflow corncob_wf {
    take:
    results_hdf
    details_hdf
    manifest_csv
    group_name
    group_col

    main:

    // Extract the counts for this grouping of genes
    // Regardless of whether a --formula is provided,
    // this will publish the table of readcounts to the output
    extractCounts(
        results_hdf,
        details_hdf,
        group_name,
        group_col
    )

    // If a formula has also been provided
    if ( params.formula ) {

        // Set up a channel with the strings of the formula(s) provided
        formula_ch = Channel.of(
            params.formula.split(",")
        )

        splitCorncob(
            extractCounts.out[0]
        )

        runCorncob(
            splitCorncob.out.flatten(),
            manifest_csv,
            group_name,
            formula_ch
        )

        joinCorncob(
            runCorncob.out.toSortedList(),
            group_name
        )

        addCorncobResults(
            extractCounts.out[2],
            joinCorncob.out,
            group_name
        )

        output_hdf = addCorncobResults.out

    } else {

        output_hdf = extractCounts.out[2]

    }

    emit:
    output_hdf
}

process mockData {
    tag "Simulate dataset for validation"
    container "${container__pandas}"
    label 'io_limited'

    input:
    path manifest_csv

    output:
    path "random.counts.csv.gz"

"""
#!/usr/bin/env python3

import numpy as np
import pandas as pd

# Get the user-provided manifest
manifest = pd.read_csv("${manifest_csv}")

# Make sure that 'specimen' is in the header
assert 'specimen' in manifest.columns.values, "Must provide a column named 'specimen'"

# Get the list of specimen names
specimen_names = list(set(manifest['specimen'].tolist()))

# Make a new DataFrame with random numbers
df = pd.DataFrame(
    np.random.randint(
        1000, high=2000, size=[len(specimen_names), 5], dtype=int
    ),
    index=specimen_names,
    columns=[
        "CAG-%d" % i
        for i in range(5)
    ]
)

# Add a total column
df["total"] = df.sum(axis=1)

# Write out to a file
df.reset_index(
).rename(
    columns=dict([("index", "specimen")])
).to_csv(
    "random.counts.csv.gz",
    index=None,
    compression="gzip"
)

"""

}

process validateFormula {
    tag "Validate user-provided formula"
    container "${container__pandas}"
    label 'io_limited'

    input:
    path corncob_output_csv
    path manifest_csv

    output:
    path "${manifest_csv}"

"""
#!/usr/bin/env python3

import pandas as pd

# Open up the corncob results for the simulated random data
df = pd.read_csv("${corncob_output_csv}")

# Make sure that we have results for every CAG
assert set(df["CAG"].tolist()) == set(["CAG-%d" % i for i in range(5)])

# Check to see if any CAGs returned 'failed'
if "failed" in df["type"].values:
    n_failed_cags = df.query("type == 'failed'")["CAG"].unique().shape[0]
    n_total_cags = df["CAG"].unique().shape[0]
    msg = "%d / %d CAGs failed processing with this formula" % (n_failed_cags, n_total_cags)
    assert False, msg

# Make sure that every CAG has each expected row
print("Making sure that every CAG has results for p_value, std_error, and estimate")
for cag_id, cag_df in df.groupby("CAG"):

    msg = "%s: Found %s" % (cag_id, ", t".join(cag_df["type"].tolist()))
    assert set(cag_df["type"].tolist()) == set(["p_value", "std_error", "estimate"]), msg
"""

}


// Extract a table with the number of reads per gene group
// This can be used to extract the counts for CAGs, tax IDs, or any
// other grouping of genes
// Corncob takes the absolute number of reads from each sample into account
// and so it needs to have access to those integer values
process extractCounts {
    container "${container__pandas}"
    label 'mem_veryhigh'
    publishDir "${params.output_folder}/abund/", mode: "copy"
    
    input:
    file results_hdf
    file details_hdf
    val group_name
    val group_col

    output:
    file "${group_name}.readcounts.csv.gz"
    file "${group_name}.abundance.csv.gz"
    file "${results_hdf}"


"""
#!/usr/bin/env python3

from collections import defaultdict
import gzip
import json
import logging
import os
import pandas as pd

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [Extract Counts ${group_name}] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Read in the group assignments for each gene

logging.info("Reading in gene groupings from ${results_hdf}")
logging.info("Reading in ${group_name} assignments")
gene_group_df = pd.read_hdf("${results_hdf}", "/annot/gene/all")

# Make sure that the expected columns are present

msg = "Expected column 'gene' to be present"
assert "gene" in gene_group_df.columns.values, msg

msg = "Expected column '${group_col}' to be present"
assert "${group_col}" in gene_group_df.columns.values, msg

# Format the gene groups as a Series (of integers)
gene_groups = gene_group_df.set_index(
    "gene"
)[
    "${group_col}"
].dropna(    
).apply(
    int
)

# Make an object to hold the number of reads per group, per sample
group_counts = defaultdict(lambda: defaultdict(int))

# Make an object to hold the depth-based relative abundance per group, per sample
group_abund = defaultdict(lambda: defaultdict(int))

# Read in each of the FAMLI output objects from the details HDF

logging.info("Reading in gene abundances from ${details_hdf}")

# Open the HDF store
with pd.HDFStore("${details_hdf}", "r") as store:

    # Iterate over all keys
    for hdf_key in store.keys():

        # If the key holds detailed abundances
        if hdf_key.startswith("/abund/gene/long/"):

            # Parse the specimen name
            specimen = hdf_key[len('/abund/gene/long/'):]

            logging.info(
                "Reading in %s" % specimen
            )

            # Read the complete table of abundances
            df = pd.read_hdf(store, hdf_key).set_index("id")

            # Compute the relative depth
            df = df.assign(rel_depth = df['depth'] / df['depth'].sum())

            # Iterate over every gene
            for gene_id, r in df.iterrows():

                # Add the number of reads to the count for the group
                # If the gene does not have any group assigned, then
                # the count will go to "UNASSIGNED", to preserve
                # the appropriate total read counts per specimen

                # First key is specimen
                group_counts[
                    specimen
                ][  # Second key is gene group
                    gene_groups.get(
                        gene_id,
                        "UNASSIGNED"
                    )
                ] += r['nreads']

                # Saving the relative abundance based on depth
                group_abund[
                    specimen
                ][  # Second key is gene group
                    gene_groups.get(
                        gene_id,
                        "UNASSIGNED"
                    )
                ] += r['rel_depth']

# Function to save a dict of dicts as a CSV, and to HDF
def save_data(
    group_dict,
    file_path,
    store,
    hdf_key,
    data_type
):

    # Format as a DataFrame
    df = pd.DataFrame(
        group_dict
    ).fillna(
        0
    ).applymap(
        data_type
    )

    # Transform so that specimens are rows
    df = df.T

    # Reset the index and add a name indicating that the rows 
    # correspond to specimens from the manifest
    df = df.reset_index(
    ).rename(
        columns=dict([("index", "specimen")])
    )

    # Save to a file
    logging.info("Writing to %s" % file_path)
    df.to_csv(
        file_path,
        compression="gzip",
        index=None
    )

    # Save to HDF
    df.to_hdf(store, hdf_key, format='fixed')

# Open the results HDF for writing
with pd.HDFStore("${results_hdf}", 'a') as store:

    logging.info("Making a DataFrame of ${group_name} counts")
    save_data(
        group_counts,
        "${group_name}.readcounts.csv.gz",
        store,
        "/counts/${group_name}/wide",
        int
    )

    save_data(
        group_abund,
        "${group_name}.abundance.csv.gz",
        store,
        "/abund/${group_name}/wide",
        float
    )


logging.info("Done")
"""
}


process splitCorncob {
    container "${container__pandas}"
    label "mem_veryhigh"
    errorStrategy "retry"

    input:
    file readcounts_csv_gz

    output:
    file "readcounts.*.csv.gz"

"""#!/usr/bin/env python3

import pandas as pd

# Read in the entire set of readcounts
print("Reading in ${readcounts_csv_gz}")
df = pd.read_csv("${readcounts_csv_gz}")
print("Read in %d rows and %d columns" % (df.shape[0], df.shape[1]))

# Add in a 'total' column
df = df.assign(
    total = df.drop(columns='specimen').sum(axis=1)
)

# Write out shards of the data
for shard_ix in range(${params.corncob_batches}):

    # Get the list of CAGs to write out
    cag_list = [
        n
        for ix, n in enumerate(df.columns.values)
        if ix % ${params.corncob_batches} == shard_ix or n in ["specimen", "total"]
    ]

    # Skip if there are too few CAGs for this shard
    if len(cag_list) <= 2:
        print("Skipping shard %d, too few CAGs to populate" % shard_ix)
        continue

    # Write out the shard
    print("Writing out %d columns to shard %d" % (len(cag_list), shard_ix))

    df.reindex(columns=cag_list).to_csv("readcounts.%d.csv.gz" % shard_ix, index=None)
    print("Done with shard %d" % shard_ix)

print("Done with all shards")
"""
}

// Extract the counts table from the results HDF5
process runCorncob {
    tag "Perform statistical analysis"
    container "quay.io/fhcrc-microbiome/corncob"
    label "mem_medium"
    errorStrategy "retry"
    
    input:
    file readcounts_csv_gz
    file metadata_csv
    val group_name
    each formula

    output:
    file "corncob.results.csv"


    """
#!/usr/bin/env Rscript

# Get the arguments passed in by the user

library(tidyverse)
library(corncob)
library(parallel)

## By default, use 10% of the available memory to read in data
connectionSize = 100000 * ${task.memory.toMega()}
print("Using VROOM_CONNECTION_SIZE =")
print(connectionSize)
Sys.setenv("VROOM_CONNECTION_SIZE" = format(connectionSize, scientific=F))

numCores = ${task.cpus}

##  READCOUNTS CSV should have columns `specimen` (first col) and `total` (last column).
##  METADATA CSV should have columns `specimen` (which matches up with `specimen` from
##         the recounts file), and additional columns with covariates matching `formula`

##  corncob analysis (coefficients and p-values) are written to OUTPUT CSV on completion

print("Reading in ${metadata_csv}")
metadata <- vroom::vroom("${metadata_csv}", delim=",")

print("Removing columns which are not in the formula")
for(column_name in names(metadata)){
    if(column_name == "specimen" || grepl(column_name, "${formula}", fixed=TRUE) ){
        print(paste("Keeping column", column_name))
    } else {
        print(paste("Removing column", column_name))
        metadata <- metadata %>% select(-column_name)
    }
}
metadata <- metadata %>% unique %>% drop_na
print("Filtered and deduplicated manifest:")
print(metadata)

print("Reading in ${readcounts_csv_gz}")
counts <- vroom::vroom("${readcounts_csv_gz}", delim=",")
total_counts <- counts[,c("specimen", "total")]

print("Adding total counts to manifest")
print(head(total_counts))

print("Merging total counts with metadata")
total_and_meta <- metadata %>% 
  left_join(total_counts, by = c("specimen" = "specimen"))

#### Run the analysis for every individual ${group_name} (in this shard)
print(sprintf("Starting to process %s columns (${group_name})", length(c(2:(dim(counts)[2] - 1)))))
corn_tib <- do.call(rbind, mclapply(
    c(2:(dim(counts)[2] - 1)),
    function(i){
        try_bbdml <- try(
            counts[,c(1, i)] %>%
            rename(W = 2) %>%
            right_join(
                total_and_meta, 
                by = c("specimen" = "specimen")
            ) %>%
            corncob::bbdml(
                formula = cbind(W, total - W) ~ ${formula},
                phi.formula = ~ 1,
                data = .
            )
        )

      if (class(try_bbdml) == "bbdml") {
        return(
            summary(
                try_bbdml
            )\$coef %>%
            as_tibble %>%
            mutate("parameter" = summary(try_bbdml)\$coef %>% row.names) %>%
            rename(
                "estimate" = Estimate,
                "std_error" = `Std. Error`,
                "p_value" = `Pr(>|t|)`
            ) %>%
            select(-`t value`) %>%
            gather(key = type, ...=estimate:p_value) %>%
            mutate("${group_name}" = names(counts)[i])
        )
      } else {
          return(
              tibble(
                  "parameter" = "all",
                  "type" = "failed", 
                  "value" = NA, 
                  "${group_name}" = names(counts)[i]
              )
          )
      }   
    },
    mc.cores = numCores
  ))

print(head(corn_tib))

print("Adding a column with the formula used here")
corn_tib <- corn_tib %>% add_column(formula = "${formula}")

print(head(corn_tib))

print(sprintf("Writing out %s rows to corncob.results.csv", nrow(corn_tib)))
write_csv(corn_tib, "corncob.results.csv")
print("Done")
    """

}

// Run the breakaway algorithm on each sample
process breakaway {
    tag "Estimate richness"
    container "quay.io/fhcrc-microbiome/breakaway"
    label "io_limited"
    errorStrategy "retry"
    
    input:
    file famli_json_gz

    output:
    file "*.breakaway.json"


"""
#!/usr/bin/env Rscript

library(jsonlite)
library(breakaway)

# Read in the input data
gene_readcounts <- fromJSON("${famli_json_gz}")\$nreads

# Run breakaway
r <- breakaway(gene_readcounts, plot = FALSE, output = FALSE, answers = TRUE)

print(r)

# Make a new output object with only the data objects which are strictly needed
output <- list(
    estimate = r\$estimate,
    error = r\$error,
    interval = r\$interval,
    reasonable = r\$reasonable,
    estimand = r\$estimand
)

# Save the results to a file
output_filename <- sub(".json.gz", ".breakaway.json", "${famli_json_gz}")
write(
    toJSON(
        output,
        force = TRUE
    ),
    file = output_filename
)

"""

}

// Collect the breakaway algorithm results for all samples
process collectBreakaway {
    tag "Join richness tables"
    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
    label "io_limited"
    errorStrategy "retry"
    publishDir "${params.output_folder}stats", mode: "copy", overwrite: true
    
    input:
    file breakaway_json_list

    output:
    file "${params.output_prefix}.breakaway.csv.gz"


"""
#!/usr/bin/env python3

import json
import pandas as pd

# Get the list of files, with the samples encoded in the file names
samples = dict([
    (fp.replace(".breakaway.json", ""), fp)
    for fp in "${breakaway_json_list}".split(" ")
])

print("Reading in breakaway results for %d samples" % len(samples))

# Function to read in breakaway results
def read_breakaway(fp):
    dat = json.load(open(fp, "r"))
    return dict([
        ("estimate", dat["estimate"][0]),
        ("error", dat["error"][0]),
        ("interval_lower", dat["interval"][0]),
        ("interval_upper", dat["interval"][1]),
        ("reasonable", dat["reasonable"][0]),
        ("estimand", dat["estimand"][0])
    ])

output = pd.DataFrame(dict([
    (sample_name, read_breakaway(fp))
    for sample_name, fp in samples.items()
])).T.reset_index(
).rename(columns=dict([("index", "specimen")]))

output.to_csv("${params.output_prefix}.breakaway.csv.gz", index=None)

"""

}

// Join together a set of corncob results CSVs
process joinCorncob {
    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
    label "io_limited"
    errorStrategy "retry"
    publishDir "${params.output_folder}/stats/", mode: "copy"
    
    input:
    file "corncob.results.*.csv"
    val group_name

    output:
    file "corncob.results.csv"


"""
#!/usr/bin/env python3

import os
import pandas as pd

# Get the list of files to join
fp_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("corncob.results.") and fp.endswith(".csv")
]

print("Reading in corncob results for %d formula(s)" % len(fp_list))

df = pd.concat([
    pd.read_csv(fp)
    for fp in fp_list
]).query(  # Filter out the unassigned groups
    "${group_name} != 'UNASSIGNED'"
)

print("Writing out to corncob.results.csv")
df.to_csv("corncob.results.csv", index=None)
print("Done")
"""

}


process addCorncobResults{
    tag "Add statistical analysis to HDF"
    container "${container__pandas}"
    label 'mem_veryhigh'
    errorStrategy 'finish'

    input:
        path results_hdf
        path corncob_csv
        val group_name

    output:
        path "${results_hdf}"

"""
#!/bin/bash

add_corncob_results.py "${results_hdf}" "${corncob_csv}" "${params.fdr_method}" "${group_name}"

"""

}