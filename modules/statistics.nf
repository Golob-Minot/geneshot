// // Processes to perform statistical analysis

// Container versions
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

// Workflow to test out the provided formula and manifest CSV in a dry run
// This is intended to catch the case early on where the provided formula
// is not formatted correctly, or does not match the manifest (sample sheet)
workflow validation_wf {
    take: 
        manifest_csv
        formula_ch
    main: 

    mockData(
        manifest_csv
    )

    runCorncob(
        mockData.out,
        manifest_csv,
        formula_ch
    )

    joinCorncob(
        runCorncob.out.toSortedList()
    )

    validateFormula(
        joinCorncob.out,
        manifest_csv
    )

    emit:
    validateFormula.out

}

// Workflow to run corncob on the actual data
workflow corncob_wf {
    take:
    famli_json_list
    cag_csv
    manifest_csv
    formula_ch

    main:

    extractCounts(
        famli_json_list,
        cag_csv
    )

    runCorncob(
        extractCounts.out,
        manifest_csv,
        formula_ch
    )

    joinCorncob(
        runCorncob.out.toSortedList()
    )

    emit:
    joinCorncob.out
}

process mockData {
    tag "Simulate dataset for validation"
    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
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
    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
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

// Extract a CAG-level counts table from the FAMLI JSON outputs
// Corncob takes the absolute number of reads from each sample into account
// and so it needs to have access to those integer values
process extractCounts {
    tag "Make CAG ~ sample read-count matrix"
    container "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
    label 'mem_veryhigh'
    publishDir "${params.output_folder}/abund/", mode: "copy"
    
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
    tag "Perform statistical analysis"
    container "quay.io/fhcrc-microbiome/corncob"
    label "mem_veryhigh"
    errorStrategy "retry"
    
    input:
    file readcounts_csv_gz
    file metadata_csv
    val formula

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
            mutate("CAG" = names(counts)[i])
        )
      } else {
          return(
              tibble(
                  "parameter" = "all",
                  "type" = "failed", 
                  "value" = NA, 
                  "CAG" = names(counts)[i]
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
])

print("Writing out to corncob.results.csv")
df.to_csv("corncob.results.csv", index=None)
print("Done")
"""

}


// Run meta-analysis on corncob results grouped by annotation label
// Each input will have the results from a single type of annotation:
// species, genus, family, or eggNOG_desc
// The column `label` gives the label, and there is a row for every 
// CAG which has at least one gene with that label
// Other columns include CAG, parameter, estimate, and std_error
process runBetta {
    container "quay.io/fhcrc-microbiome/breakaway:latest"
    label "mem_medium"
    errorStrategy "retry"
    
    input:
    path labelled_corncob_csv

    output:
    file "${labelled_corncob_csv}.betta.csv.gz"


"""
#!/usr/bin/env Rscript
library(tidyverse)
library(magrittr)
library(reshape2)
library(breakaway)

# Use vroom to read in the table
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 20)

# Read in all of the data for a single covariate
print("Reading in $labelled_corncob_csv")
df <- vroom::vroom("$labelled_corncob_csv", delim=",")
df <- as.tibble(df)

print(head(df))

# Make a function which will trim off the trailing digit
options(scipen=999)
num_decimal_places <- function(v){nchar(strsplit(as.character(v), "\\\\.")[[1]][2])}
trim_trailing <- function(v){round(v, num_decimal_places(v) - 1)}

# Make a function to run betta while tolerating faults
fault_tolerant_betta <- function(df, f){
    if(nrow(df) == 1){
        return(
            data.frame(
                estimate=df[1,"estimate"],
                std_error=df[1,"std_error"],
                p_value=df[1,"p_value"]
            )
        )
    }
    chats <- df\$estimate
    ses <- df\$std_error
    r <- NULL
    for(ix in c(1:10)){
        r <- tryCatch({
            breakaway::betta(chats=chats, ses=ses)\$table
            },
            error=function(cond) {
                print("We have encountered an error:")
                print(cond)
                print("The input data which caused the error was:")
                print(chats)
                print(ses)
                return(NULL)
            })
        if(!is.null(r)){
            return(
                data.frame(
                    estimate=r[1,"Estimates"],
                    std_error=r[1,"Standard Errors"],
                    p_value=r[1,"p-values"]
                )
            )
        } else {
            print("Trimming down the input data")
            chats <- trim_trailing(chats)
            ses <- trim_trailing(ses)
            print(chats)
            print(ses)
        }
    }
}

# If there is a single dummy row, skip the entire process
if(nrow(df) == 1){

    write.table(df, file=gzfile("${labelled_corncob_csv}.betta.csv.gz"), sep=",", row.names=FALSE)

} else{

    # Perform meta-analysis combining the results for each label, and each parameter
    results <- df %>% group_by(annotation, label, parameter) %>% group_modify(fault_tolerant_betta)

    # Write out to a CSV
    write.table(results, file=gzfile("${labelled_corncob_csv}.betta.csv.gz"), sep=",", row.names=FALSE)

}

"""

}


process addBetta{
    tag "Add meta-analysis to HDF"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
        path results_hdf
        path betta_csv_list

    output:
        path "${results_hdf}"

"""
#!/usr/bin/env python3

import os
import pandas as pd
from statsmodels.stats.multitest import multipletests

betta_csv_list = "${betta_csv_list}".split(" ")

for betta_csv in betta_csv_list:
    if len(betta_csv) > 1:
        assert os.path.exists(betta_csv)

# Read in from the flat file
df = pd.concat([
    pd.read_csv(betta_csv)
    for betta_csv in betta_csv_list
    if len(betta_csv) > 1
])

print("Read in {:,} lines from {}".format(
    df.shape[0],
    betta_csv
))

# If there are real results (not just a dummy file), write to HDF5
if df.shape[0] > 1:

    # Add the q-values
    df = df.assign(
        q_value = multipletests(
            df["p_value"],
            0.2,
            "${params.fdr_method}"
        )[1]
    )

    # Open a connection to the HDF5
    with pd.HDFStore("${results_hdf}", "a") as store:

        # Write to HDF5
        key = "/stats/enrichment/betta"
        print("Writing to %s" % key)
        
        # Write to HDF
        df.to_hdf(store, key)

        print("Closing store")

    print("Done")

else:
    print("No betta results found -- returning unopened results HDF")

"""

}
