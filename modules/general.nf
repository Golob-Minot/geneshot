
process combineReads {
    container 'quay.io/fhcrc-microbiome/fastatools@sha256:171884eba1fb5561a8a42cf23a5ffd9857e6e34e1331d9e5dca220ada44555a7'
    label = 'io_limited'
    errorStrategy 'retry'
    maxRetries 10

    // If the user sets --preprocess_output, write out the combined reads to that folder
    publishDir path: "${params.output_folder}qc/", enabled: params.savereads, mode: "copy"

    input:
    tuple val(sample), file("R1.*.fastq.gz"), file("R2.*.fastq.gz")
    
    output:
    tuple val(sample), file("${sample}.R1.fastq.gz"), file("${sample}.R2.fastq.gz")

"""
set -e

ls -lah *

combine_fastq_pairs.py \
-1 R1*fastq.gz \
-2 R2*fastq.gz \
--normalize-ids \
-o1 "${sample}.R1.fastq.gz" \
-o2 "${sample}.R2.fastq.gz"

(( \$(gunzip -c "${sample}.R1.fastq.gz" | wc -l) > 1 ))
(( \$(gunzip -c "${sample}.R2.fastq.gz" | wc -l) > 1 ))

"""

}

process outputManifest {
    container "ubuntu:18.04"

    publishDir path: "${params.output_folder}qc/", enabled: params.savereads, mode: "copy"

    input:
        val manifestStr
    
    output:
        file 'manifest.qc.csv'

    """
        echo "${manifestStr}" > manifest.qc.csv
    """
}

// Workflow to publish a set of reads to a folder, along with a manifest
workflow writeManifest {
    get:
        reads_ch

    main:
        // Make a manifest for the files in reads_ch
        // Output the final reads and manifest

        
        manifestStr = reads_ch.reduce(
            'specimen,R1,R2\n'
        ){ csvStr, row ->
            return  csvStr += "${row[0]},${params.output_folder}qc/${row[1].name},${params.output_folder}qc/${row[2].name}\n";
        }

        // Write the manifest CSV to a file
        outputManifest(manifestStr)
        
}

// Count the number of input reads for a single sample
process countReads {
    container "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.1A"
    cpus 1
    memory "4 GB"
    errorStrategy "retry"

    input:
    tuple sample_name, file(fastq)

    output:
    file "${sample_name}.countReads.csv"

"""
set -e

[[ -s ${fastq} ]]

n=\$(gunzip -c "${fastq}" | awk 'NR % 4 == 1' | wc -l)
echo "${sample_name},\$n" > "${sample_name}.countReads.csv"
"""
}


// Make a single file which summarizes the number of reads across all samples
// This is only run after all of the samples are done processing through the
// 'total_counts' channel, which is transformed by the .collect() command into
// a single list containing all of the data from all samples.
process countReadsSummary {
    container "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.1A"
    // The output from this process will be copied to the --output_folder specified by the user
    publishDir "${params.output_folder}/qc/", mode: 'copy'
    errorStrategy "retry"

    input:
    // Because the input channel has been collected into a single list, this process will only be run once
    file readcount_csv_list

    output:
    file "readcounts.csv"


"""
set -e

echo name,n_reads > readcounts.csv
cat ${readcount_csv_list} >> readcounts.csv
"""
}


// Process which will concatenate a set of files
process concatenateFiles {
    container "ubuntu:18.04"
    label "mem_medium"
    errorStrategy "retry"
    
    input:
    file "__INPUT*"
    val output_name

    output:
    file "${output_name}"

"""
# Break on any errors
set -e

cat __INPUT* > ${output_name}
"""
}

process collectAbundances{
    container "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
        path cag_csv
        path gene_abund_feather
        path cag_abund_feather
        path famli_json_list
        path readcount_csv

    output:
        path "results.hdf5"

"""
#!/usr/bin/env python3

import gzip
import json
import pandas as pd

cag_csv = "${cag_csv}"
gene_abund_feather = "${gene_abund_feather}"
cag_abund_feather = "${cag_abund_feather}"
famli_json_list = "${famli_json_list}".split(" ")
readcount_csv = "${readcount_csv}"

# Function to read in the FAMLI output
def read_famli_json(fp, suffix=".json.gz"):

    # Get the sample name from the file name
    assert fp.endswith(suffix)
    sample_name = fp[:-len(suffix)]

    return pd.DataFrame(
        json.load(
            gzip.open(
                fp, "rt"
            )
        )
    ).assign(
        specimen=sample_name
    )


# Open a connection to the HDF5
with pd.HDFStore("results.hdf5", "w") as store:

    # Read in the complete set of FAMLI results
    famli_df = pd.concat([
        read_famli_json(fp)
        for fp in famli_json_list
    ])

    # Write to HDF5
    famli_df.to_hdf(store, "/abund/gene/long")
    
    # Read in the summary of the number of reads across all samples
    readcount_df = pd.read_csv(readcount_csv)

    # Write to HDF5
    readcount_df.to_hdf(store, "/summary/readcount")
    

    # Read in the table with the CAG-level abundances across all samples
    cag_abund_df = pd.read_feather(
        cag_abund_feather
    )
    print(
        "Read in abundances for %d CAGs across %d samples" %
        (cag_abund_df.shape[0], cag_abund_df.shape[1] - 1)
    )


    # Write to HDF5
    cag_abund_df.to_hdf(store, "/abund/cag/wide")

    # Read in the table with the gene-level abundances across all samples
    gene_abund_df = pd.read_feather(
        gene_abund_feather
    )
    print(
        "Read in abundances for %d genes across %d samples" %
        (gene_abund_df.shape[0], gene_abund_df.shape[1] - 1)
    )

    # Write to HDF5
    gene_abund_df.to_hdf(store, "/abund/gene/wide")

    # Read in the table describing which genes are grouped into which CAGs
    cag_df = pd.read_csv(cag_csv)

    print(
        "Read in CAG assignments for %d genes across %d CAGs" % 
        (cag_df.shape[0], cag_df["CAG"].unique().shape[0])
    )

    # Write to HDF5
    cag_df.to_hdf(store, "/annot/gene/cag")

"""

}

process addGeneAssembly{
    container "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
        path results_hdf
        path allele_gene_tsv
        path allele_assembly_csv

    output:
        path "${results_hdf}"

"""
#!/usr/bin/env python3

import pandas as pd

# Read in the summary of allele assembly
allele_assembly = pd.read_csv("${allele_assembly_csv}")

print(
    "Read in %d gene assemblies over %d specimens" % 
    (allele_assembly.shape[0], allele_assembly["specimen"].unique().shape[0])
)

# Read in the listing of which alleles were grouped into which genes
allele_gene = pd.read_csv(
    "${allele_gene_tsv}", 
    sep="\\t", 
    header=None
).loc[
    :, :1
].rename(
    columns=dict([
        (0, "gene"),
        (1, "allele")
    ])
)
print(
    "Read in a table grouping %d alleles into %d genes" % 
    (
        allele_gene["allele"].unique().shape[0],
        allele_gene["gene"].unique().shape[0]
    )
)

# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Write assembly summary to HDF5
    allele_assembly.to_hdf(store, "/abund/allele/assembly")

    # Write gene <-> allele table to HDF5
    allele_gene.to_hdf(store, "/annot/allele/gene")

"""

}

process addCorncobResults{
    container "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
        path results_hdf
        path corncob_csv

    output:
        path "${results_hdf}"

"""
#!/usr/bin/env python3

import pandas as pd

# Read in the corncob results
corncob_df = pd.read_csv("${corncob_csv}")

print(
    "Read in corncob results for %d CAGs" % 
    corncob_df["CAG"].unique().shape[0]
)


# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Write corncob results to HDF5
    corncob_df.to_hdf(store, "/stats/cag/corncob")

"""

}

// Repack an HDF5 file
process repackHDF {

    container "quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed"
    label "mem_veryhigh"
    errorStrategy "retry"
    
    input:
    file output_hdf5
        
    output:
    file "${output_hdf5}"

    """
#!/bin/bash

set -e

[ -s ${output_hdf5} ]

h5repack -f GZIP=5 ${output_hdf5} TEMP && mv TEMP ${output_hdf5}
    """
}