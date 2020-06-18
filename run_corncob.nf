#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

  This utility runs the corncob algorithm with a specified formula on an
  existing geneshot output file object.
*/

// Using DSL-2
nextflow.preview.dsl=2

// Parameters
params.input_hdf = false
params.input_folder = false
params.output_folder = false
params.output_prefix = false
params.formula = false
params.fdr_method = "fdr_bh"
params.help = false

// Docker containers
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Run corncob with a specified formula on a geneshot output file.
    This utility will create a new geneshot output file in which
    the only difference is the corncob output results stored within.

    Usage:

    nextflow run Golob-Minot/geneshot/run_corncob.nf <ARGUMENTS>
    
    Options:
      --input_hdf           Geneshot results HDF5 file to be used for analysis
      --input_folder        Folder containing geneshot output to be used for analysis (must contain the subfolder "abund/")
      --output_folder       Folder to place output tile
      --output_prefix       Text used as a prefix for summary HDF5 output files (the .corncob.hdf5 suffix will be attached)
      -w                    Working directory. Defaults to `./work`
      
    For Statistical Analysis:
      --formula             Optional formula used to estimate associations with CAG relative abundance
      --fdr_method          FDR method used to calculate q-values for associations (default: 'fdr_bh')

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.input_hdf == false || params.input_folder == false || params.output_folder == false || params.output_prefix == false || params.formula == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Import the process used for statistical analysis
include runCorncob from './modules/statistics' params(
    formula: params.formula
)

// Import the process used to join the results for multiple formulae
include joinCorncob from './modules/statistics' params(
    output_folder: params.output_folder
)

// Import the process used to aggregate readcounts across all samples
// This process is only run if `abund/CAG.readcounts.csv.gz` can not be found yet
// After the process completes, a new `abund/CAG.readcounts.csv.gz` will be added to that input folder
include extractCounts from './modules/statistics' params(
    output_folder: params.input_folder
)

// Import the processes needed to run meta-analysis of corncob results by annotation
include runBetta from './modules/statistics'
include addBetta from './modules/statistics'

// Import the process used to add corncob results to the output
include repackHDF from './modules/general' params(
    output_folder: params.output_folder
)
include addCorncobResults from './modules/general' params(
    fdr_method: params.fdr_method
)

// Process to update the formula listed in the summary table
// Also extract the manifest to reduce the number of times
// the file is opened, as well as the table listing the genes
// which make up each CAG.
process updateFormula{
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
        path results_hdf

    output:
        path "*corncob.hdf5"
        path "manifest.csv"
        path "CAG.assignments.csv.gz"

"""
#!/usr/bin/env python3

import pandas as pd
import shutil

# Come up with a new name for the output file
formula_name = "${params.formula}".replace(
    " ", "_"
).replace(
    "+", "_"
).replace(
    ",", "_"
).replace(
    "*", "_"
).replace(
    ":", "_"
).replace(
    "__", "_"
).replace(
    "__", "_"
)
new_file_name = "${params.output_prefix}.%s.corncob.hdf5" % formula_name

# Copy the input file to a new output file object
print("Copying %s to %s" % ("${results_hdf}", new_file_name))
shutil.copyfile("${results_hdf}", new_file_name)

# Open a connection to the store
with pd.HDFStore(new_file_name, "a") as store:
    
    # Read in the existing summary table
    df = pd.read_hdf(store, "/summary/experiment")

    print("Existing summary:")
    print(df)

    # Add the formula
    df = pd.concat([
        df.query("variable != 'formula'"),
        pd.DataFrame([dict(variable='formula', value='${params.formula}')])
    ])

    print("")
    print("New summary:")
    print(df)

    # Write back to the store
    df.to_hdf(store, "/summary/experiment")

    # Extract the manifest
    print("Extracting the manifest")
    pd.read_hdf(store, "/manifest").to_csv("manifest.csv")

    # Extract the table of genes making up each CAG
    print("Extracting the gene~CAG table")
    pd.read_hdf(
        store, 
        "/annot/gene/all",
        columns=["gene", "CAG"]
    ).to_csv(
        "CAG.assignments.csv.gz",
        index=None
    )

print("Done")
"""

}

workflow {

    // Make sure we can find the input file
    if(file(params.input_hdf).isEmpty()){
        log.info"""Cannot find input file ${params.input_hdf}""".stripIndent()
        exit 0
    }

    // First update the formula listed in the geneshot results file,
    // a process which will also extract the manifest CSV
    updateFormula(
        file(params.input_hdf)
    )

    // Check to see if `abund/CAG.readcounts.csv.gz` is present in the input folder
    // For context, this file is only created by geneshot for runs which originally
    // included a `--formula`. 
    // If this file cannot be found, then we will run a process to create that file
    // and also to publish that file back to the input folder.
    abund_folder = "${params.input_folder.replaceAll(/\/$/, "")}/abund"
    if(file("${abund_folder}/CAG.readcounts.csv.gz").isEmpty()){
        extractCounts(
            Channel.fromPath(
                "${abund_folder}/details/*json.gz"
            ).toSortedList(),
            updateFormula.out[2]
        )
        input_csv = extractCounts.out
    } else {
        input_csv = file("${abund_folder}/CAG.readcounts.csv.gz")
    }

    // Set up a channel with the strings of the formula(s) provided
    formula_ch = Channel.of(
        params.formula.split(",")
    )

    // Now run corncob on the extracted manifest, as well as the gene counts table
    runCorncob(
        input_csv,
        updateFormula.out[1],
        formula_ch
    )

    joinCorncob(
        runCorncob.out.toSortedList()
    )

    // Add those results to the output file
    addCorncobResults(
        updateFormula.out[0],
        joinCorncob.out
    )

    runBetta(
        addCorncobResults.out[1].flatten()
    )

    addBetta(
        runBetta.out.toSortedList(),
        addCorncobResults.out[0]
    )

    // Repack the HDF
    repackHDF(
        addBetta.out[0]
    )

}