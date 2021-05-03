#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
  are present in a microbial community.

  This utility replaces the /manifest table in an existing geneshot results output file.
*/

// Using DSL-2
nextflow.enable.dsl=2

// Parameters
params.input_hdf = false
params.output_folder = false
params.output_prefix = false
params.manifest = false
params.help = false

// Docker containers
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Replace the /manifest in an existing geneshot output file.

    Usage:

    nextflow run Golob-Minot/geneshot/update_manifest.nf <ARGUMENTS>
    
    Options:
      --manifest            Manifest (CSV) to be used in the updated geneshot output file
      --input_hdf           Geneshot results HDF5 file to be used for analysis
      --output_folder       Folder to place output tile
      --output_prefix       Text used as a prefix for summary HDF5 output files (the .corncob.hdf5 suffix will be attached)
      -w                    Working directory. Defaults to `./work`
      
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.input_hdf == false || params.manifest == false || params.output_folder == false || params.output_prefix == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Import the process used to repack the output HDF5 file
include repackHDF from './modules/general'

// Process to update the formula listed in the summary table
// Also extract the manifest to reduce the number of times
// the file is opened
process replaceManifest {
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'finish'

    input:
        path "input.hdf5"
        path manifest_csv

    output:
        path "${params.output_prefix}.results.hdf5"

"""
#!/usr/bin/env python3

import pandas as pd
import shutil
import pickle
pickle.HIGHEST_PROTOCOL = 4

# Read in the new manifest
new_manifest = pd.read_csv("${manifest_csv}", sep=",")

# Make sure that the new manifest has the three required columns
for k in ["specimen", "R1", "R2"]:
    assert k in new_manifest.columns.values, "New manifest must contain a %s column" % k

# Remove the read path columns, if present
for k in ["R1", "R2", "I1", "I2"]:
    if k in new_manifest.columns.values:
        new_manifest = new_manifest.drop(columns=k)

# Assign the new name for the output file
new_file_name = "${params.output_prefix}.results.hdf5"

# Copy the input file to a new output file object
print("Copying input.hdf5 to %s" % new_file_name)
shutil.copyfile("input.hdf5", new_file_name)

# Open a connection to the store
with pd.HDFStore(new_file_name, "a") as store:
    
    # Read in the existing manifest
    old_manifest = pd.read_hdf(store, "/manifest")

    # Make sure that the new and old manifest have the same set of specimens
    old_specimen_set = set(old_manifest["specimen"].tolist())
    new_specimen_set = set(new_manifest["specimen"].tolist())

    for s, msg in [
        (old_specimen_set - new_specimen_set, "Specimens missing from old specimen list"),
        (new_specimen_set - old_specimen_set, "Specimens missing from new specimen list"),
    ]:
        assert len(s) == 0, "%s: %s" % (msg, ", ".join(list(s)))

    # Write new manifest to the store
    new_manifest.to_hdf(store, "/manifest")

print("Done")
"""

}

workflow {

    // Update the manifest
    replaceManifest(
        file(params.input_hdf),
        file(params.manifest)
    )

    // Repack the HDF
    repackHDF(
        replaceManifest.out
    )

    // Publish output files
    publish:
        repackHDF.out to: "${params.output_folder}", mode: "copy", overwrite: true

}