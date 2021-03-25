#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

  This utility runs the corncob algorithm with a specified formula on an
  existing geneshot output file object.
*/

// Using DSL-2
nextflow.enable.dsl=2

// Parameters
params.results_hdf = false
params.details_hdf = false
params.output_folder = false
params.output_prefix = false
params.formula = false
params.fdr_method = "fdr_bh"
params.corncob_batches = 10
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
      --results_hdf         Geneshot results HDF5 file to be used for analysis
      --details_hdf         Geneshot details HDF5 file to be used for analysis
      --output_folder       Folder to place output file
      --output_prefix       Text used as a prefix for summary HDF5 output files (the .corncob.hdf5 suffix will be attached)
      -w                    Working directory. Defaults to `./work`
      
    For Statistical Analysis:
      --formula             Optional formula used to estimate associations with CAG relative abundance
      --fdr_method          FDR method used to calculate q-values for associations (default: 'fdr_bh')
      --corncob_batches     Number of parallel processes to use processing each formula

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.results_hdf == false || params.details_hdf == false || params.output_folder == false || params.output_prefix == false || params.formula == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Import the workflows used for statistical analysis
// Use separate workflows for corncob, each
// to analyze either the CAGs, tax IDs, or eggNOG groups
include { 
    validation_wf; 
    corncob_wf as cag_corncob_wf;
    corncob_wf as tax_corncob_wf;
    corncob_wf as func_corncob_wf;
} from './modules/statistics' params(
    output_folder: params.output_folder,
    formula: params.formula,
    corncob_batches: params.corncob_batches,
    fdr_method: params.fdr_method,
    output_prefix: params.output_prefix,    
)



// Import the process used to compress the output HDF
include {
    repackHDF
 } from './modules/general' params(
    output_folder: params.output_folder,
)

// Process to update the formula listed in the summary table
// Also extract the manifest to reduce the number of times
// the file is opened
process updateFormula{
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'finish'

    input:
        path results_hdf

    output:
        path "*corncob.hdf5"
        path "manifest.csv"

"""
#!/usr/bin/env python3

import pandas as pd
import shutil
import pickle
pickle.HIGHEST_PROTOCOL = 4

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

print("Done")
"""

}

workflow {

    // Make sure we can find the input file
    if(file(params.results_hdf).isEmpty()){
        log.info"""Cannot find input file ${params.results_hdf}""".stripIndent()
        exit 0
    }

    // First update the formula listed in the geneshot results file,
    // a process which will also extract the manifest CSV
    updateFormula(
        file(params.results_hdf)
    )

    // Calculate the association of individual CAGs with user-provided features
    cag_corncob_wf(
        updateFormula.out[0],
        file(params.details_hdf),
        updateFormula.out[1],
        "cag",
        "CAG"
    )

    // Repack the HDF
    repackHDF(
        cag_corncob_wf.out
    )

}