// // Processes to perform statistical analysis

workflow corncob_wf {
    get:
    famli_json_list
    cag_csv
    manifest_csv

    main:

    extractCounts(
        famli_json_list,
        cag_csv
    )

    runCorncob(
        extractCounts.out,
        manifest_csv
    )

    emit:
    runCorncob.out
}


// Extract a CAG-level counts table from the FAMLI JSON outputs
// Corncob takes the absolute number of reads from each sample into account
// and so it needs to have access to those integer values
process extractCounts {
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    label 'mem_veryhigh'
    
    input:
    file famli_json_list
    file cag_csv

    output:
    file "CAG.readcounts.csv.gz"


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
    '%(asctime)s %(levelname)-8s [extractCounts] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Read in the CAG assignment for each gene
logging.info("Reading in CAG assignments")
cags = pd.read_csv(
    "${cag_csv}"
).set_index(
    "gene"
)["CAG"]

# Make an object to hold the number of reads per CAG, per sample
cag_counts = dict()

# Read in each of the FAMLI output objects
for fp in "${famli_json_list}".split(" "):

    logging.info("Processing %s" % fp)

    # Get the sample name
    assert fp.endswith(".json.gz"), fp
    sample_name = fp[:-len(".json.gz")]

    # Make sure the file was staged correctly
    assert os.path.exists(fp), "%s not found" % fp

    # Add the counts
    cag_counts[sample_name] = defaultdict(int)

    for i in json.load(
        gzip.open(
            fp,
            "rt"
        )
    ):

        # Add the value in 'nreads' to the CAG assigned to this gene
        cag_counts[
            sample_name
        ][
            cags[
                i[
                    "id"
                ]
            ]
        ] += i[
            "nreads"
        ]

# Format as a DataFrame
logging.info("Making a DataFrame")
cag_counts = pd.DataFrame(
    cag_counts
).fillna(
    0
).applymap(
    int
)

# Transform so that samples are rows
cag_counts = cag_counts.T

# Add a "total" column
cag_counts["total"] = cag_counts.sum(axis=1)

# Reset the index and add a name indicating that the rows 
# correspond to specimens from the manifest
cag_counts = cag_counts.reset_index(
).rename(
    columns=dict([("index", "specimen")])
)

# Save to a file
logging.info("Writing to disk")
cag_counts.to_csv(
    "CAG.readcounts.csv.gz",
    compression="gzip",
    index=None
)

logging.info("Done")
"""
}

// Extract the counts table from the results HDF5
process runCorncob {
    container "quay.io/fhcrc-microbiome/corncob"
    label "mem_veryhigh"
    // errorStrategy "retry"
    
    input:
    file readcounts_csv_gz
    file metadata_csv

    output:
    file "corncob.results.csv"


    """
#!/usr/bin/env Rscript

# Get the arguments passed in by the user

library(tidyverse)
library(corncob)
library(parallel)

Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 20)

numCores = ${task.cpus}

##  READCOUNTS CSV should have columns `specimen` (first col) and `total` (last column).
##  METADATA CSV should have columns `specimen` (which matches up with `specimen` from
##         the recounts file), and additional columns with covariates matching `formula`

##  corncob analysis (coefficients and p-values) are written to OUTPUT CSV on completion

print("Reading in ${metadata_csv}")
metadata <- vroom::vroom("${metadata_csv}", delim=",")

print("Reading in ${readcounts_csv_gz}")
counts <- vroom::vroom("${readcounts_csv_gz}", delim=",")
total_counts <- counts[,c("specimen", "total")]
print("Merging total counts with metadata")
total_and_meta <- metadata %>% 
  right_join(total_counts, by = c("specimen" = "specimen"))


#### Run the analysis for every individual CAG
print(sprintf("Starting to process %s columns (CAGs)", dim(counts)[2]))
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
                formula = cbind(W, total - W) ~ ${params.formula},
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
            mutate("cag" = names(counts)[i])
        )
      } else {
          return(
              tibble(
                  "parameter" = "all",
                  "type" = "failed", 
                  "value" = NA, 
                  "cag" = names(counts)[i]
              )
          )
      }   
    },
    mc.cores = numCores
  ))

print(sprintf("Writing out %s rows to corncob.results.csv", nrow(corn_tib)))
write_csv(corn_tib, "corncob.results.csv")
    """

}
