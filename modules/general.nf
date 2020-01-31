
// Container versions
container__fastatools = "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.2"
container__ubuntu = "ubuntu:18.04"
container__experiment_collection = "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed"

process combineReads {
    tag "Join FASTQ files per-specimen"
    container "${container__fastatools}"
    label = 'mem_medium'
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
    container "${container__ubuntu}"

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
    tag "Count the number of reads per sample"
    container "${container__fastatools}"
    cpus 1
    memory "4 GB"
    errorStrategy "retry"

    input:
    tuple sample_name, file(R1), file(R2)

    output:
    file "${sample_name}.countReads.csv"

"""
set -e

[[ -s ${R1} ]]
[[ -s ${R2} ]]

n=\$(cat <(gunzip -c "${R1}") <(gunzip -c "${R2}") | awk 'NR % 4 == 1' | wc -l)
echo "${sample_name},\$n" > "${sample_name}.countReads.csv"
"""
}


// Make a single file which summarizes the number of reads across all samples
// This is only run after all of the samples are done processing through the
// 'total_counts' channel, which is transformed by the .collect() command into
// a single list containing all of the data from all samples.
process countReadsSummary {
    tag "Summarize the number of reads per sample"
    container "${container__fastatools}"
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

echo specimen,n_reads > readcounts.csv
cat ${readcount_csv_list} >> readcounts.csv
"""
}


// Process which will concatenate a set of files
process concatenateFiles {
    tag "Directly combine a group of files"
    container "${container__ubuntu}"
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
    tag "Add gene abundances to HDF"
    container "${container__experiment_collection}"
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
        path cag_csv
        path gene_abund_feather
        path cag_abund_feather
        path famli_json_list
        path readcount_csv
        path manifest_csv

    output:
        path "${params.output_prefix}.full.hdf5"

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
manifest_csv = "${manifest_csv}"

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
with pd.HDFStore("${params.output_prefix}.full.hdf5", "w") as store:

    # Read in the complete set of FAMLI results
    famli_df = pd.concat([
        read_famli_json(fp)
        for fp in famli_json_list
    ])

    # Write to HDF5
    famli_df.to_hdf(store, "/abund/gene/long")
    
    # Read in the summary of the number of reads across all samples
    readcount_df = pd.read_csv(readcount_csv)

    # Add the number of aligned reads per sample
    readcount_df["aligned_reads"] = readcount_df["specimen"].apply(
        famli_df.groupby("specimen")["nreads"].sum().get
    )

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

    # Also create the `/annot/gene/all` table, which will be added to later
    cag_df.to_hdf(store, "/annot/gene/all")

    # Write out the manifest provided by the user
    pd.read_csv("${manifest_csv}").to_hdf(store, "/manifest")

"""

}

process addGeneAssembly{
    tag "Add gene assembly data to HDF"
    container "${container__experiment_collection}"
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
    tag "Add statistical analysis to HDF"
    container "${container__experiment_collection}"
    label 'mem_medium'
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

process addEggnogResults {
    tag "Add functional predictions to HDF"
    container "${container__experiment_collection}"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
        path results_hdf
        path eggnog_csv

    output:
        path "${results_hdf}"

"""
#!/usr/bin/env python3

import pandas as pd

# Read in the eggNOG-mapper results
eggnog_df = pd.read_csv(
    "${eggnog_csv}", 
    header=3, 
    sep="\\t"
).rename(columns=dict([("#query_name", "query_name")]))

print(
    "Read in eggnog results for %d genes" % 
    eggnog_df.shape[0]
)

# Read in the existing set of gene annotations
with pd.HDFStore("${results_hdf}", "r") as store:

    # Set an index on the `gene` column
    gene_annot = pd.read_hdf(
        "${results_hdf}", 
        "/annot/gene/all"
    ).set_index("gene")

# Add in a subset of the eggnog results
gene_annot = gene_annot.assign(
    eggNOG_ortholog = eggnog_df.set_index("query_name")["seed_eggNOG_ortholog"]
).assign(
    eggNOG_tax = eggnog_df.set_index("query_name")["best_tax_level"]
).assign(
    eggNOG_desc = eggnog_df.set_index("query_name")["eggNOG free text desc."]
).reset_index()

# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Write eggnog results to HDF5
    eggnog_df.to_hdf(store, "/annot/gene/eggnog")
    
    # Write summary annotation table to HDF5
    gene_annot.to_hdf(
        store, 
        "/annot/gene/all"
    )

"""

}

process readTaxonomy {
    tag "Read the NCBI taxonomy"
    container "${container__experiment_collection}"
    label 'io_limited'
    errorStrategy 'retry'

    input:
        path taxdump_tar_gz

    output:
        path "ncbi_taxonomy.csv.gz"

"""
#!/usr/bin/env python3

import gzip
from io import BytesIO
import pandas as pd
import tarfile

tar = tarfile.open(
    "${taxdump_tar_gz}", 
    "r:gz", 
    encoding='utf-8'
)

# Read in the table of merged taxids
merged_df = pd.read_csv(
    BytesIO(
        tar.extractfile(
            "merged.dmp"
        ).read()
    ),
    sep="\\t",
    header=None
).rename(columns=dict([
    (0, "old"),
    (2, "new"),
    (6, "level")
])).reindex(
    columns=["old", "new"]
)

# Get the name of each taxon
names_df = pd.read_csv(
    BytesIO(
        tar.extractfile(
            "names.dmp"
        ).read()
    ),
    sep="\\t",
    header=None
).rename(columns=dict([
    (0, "tax_id"),
    (2, "name"),
    (6, "level")
])).query(
    "level == 'scientific name'"
).reindex(
    columns=["tax_id", "name"]
)

# Get the rank and parent of every taxon
nodes_df = pd.read_csv(
    BytesIO(
        tar.extractfile(
            "nodes.dmp"
        ).read()
    ),
    sep="\\t",
    header=None
).rename(columns=dict([
    (0, "tax_id"),
    (2, "parent"),
    (4, "rank")
])).reindex(
    columns=["tax_id", "parent", "rank"]
)

# Join the names and the nodes
tax_df = pd.concat([
    names_df.set_index("tax_id"),
    nodes_df.set_index("tax_id")
], axis=1).reset_index()

# Add in the merged taxids
tax_df = pd.concat([
    tax_df,
    tax_df.apply(
        lambda v: v.apply(merged_df.set_index("new")["old"].get) if v.name == "tax_id" else v
    ).dropna()
])

# Write to CSV
tax_df.to_csv(
    "ncbi_taxonomy.csv.gz",
    compression="gzip",
    index=None
)
"""

}


process addTaxResults {
    tag "Add taxonomic annotations to HDF"
    container "${container__experiment_collection}"
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
        path results_hdf
        path diamond_tax_csv
        path taxonomy_csv

    output:
        path "${results_hdf}"

"""
#!/usr/bin/env python3

import pandas as pd

# Read in the DIAMOND-tax results
tax_df = pd.read_csv(
    "${diamond_tax_csv}", 
    sep="\\t", 
    header=None
).rename(
    columns=dict([
        (0, "gene"), 
        (1, "tax_id"), 
        (2, "evalue")
    ])
)

print(
    "Read in taxonomic results for %d genes" % 
    tax_df.shape[0]
)

# Read in the existing set of gene annotations
with pd.HDFStore("${results_hdf}", "r") as store:

    # Set an index on the `gene` column
    gene_annot = pd.read_hdf(
        "${results_hdf}", 
        "/annot/gene/all"
    ).set_index("gene")

# Add the taxonomic annotation results
gene_annot = gene_annot.assign(
    tax_id = tax_df.set_index("gene")["tax_id"].apply(str)
).reset_index()

# Read the taxonomy CSV
taxonomy_df = pd.read_csv(
    "${taxonomy_csv}"
).applymap(
    str
)

# Add the name of the assigned tax_id
gene_annot = gene_annot.assign(
    tax_name = gene_annot["tax_id"].apply(
        taxonomy_df.set_index("tax_id")["name"].get
    )
).assign(
    tax_rank = gene_annot["tax_id"].apply(
        taxonomy_df.set_index("tax_id")["rank"].get
    )
)


# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Write taxonomic results to HDF5
    tax_df.to_hdf(store, "/annot/gene/tax")

    # Write summary gene annotation table to HDF5
    gene_annot.to_hdf(
        store, 
        "/annot/gene/all"
    )

    # Write the taxonomy table
    taxonomy_df.to_hdf(
        store,
        "/ref/taxonomy"
    )

"""

}


process makeSummaryHDF {
    tag "Make a summary HDF store"
    container "${container__experiment_collection}"
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
        path full_results_hdf

    output:
        path "${params.output_prefix}.summary.hdf5"

"""
#!/usr/bin/env python3

import pandas as pd

# Read in all of the tables needed from the full results HDF5
with pd.HDFStore("${full_results_hdf}", "r") as store:

    # Read in the manifest
    manifest_df = pd.read_hdf(store, "/manifest")

    # Table with the abundance of every CAG in every sample
    cag_abund_df = pd.read_hdf(store, "/abund/cag/wide")

    # Table with gene annotations (including CAG membership)
    gene_annot_df = pd.read_hdf(store, "/annot/gene/all")

    # Number of read pairs
    readcount_df = pd.read_hdf(store, "/summary/readcount")

    # Check to see if we have any corncob results
    if "/stats/cag/corncob" in store:

        # Read in the table of corncob results
        corncob_df = pd.read_hdf(store, "/stats/cag/corncob")
    else:
        corncob_df = None

# Now write to the new HDF store
with pd.HDFStore("${params.output_prefix}.summary.hdf5", "w") as store:

    # Write out each table, while indexing the appropriate columns
    manifest_df.to_hdf(
        store,
        "/manifest"
    )
    readcount_df.to_hdf(
        store,
        "/summary/readcount"
    )
    cag_abund_df.to_hdf(
        store, 
        "/abund/cag/wide", 
        format="table",
        data_columns=["CAG"]
    )
    gene_annot_df.to_hdf(
        store, 
        "/annot/gene/all", 
        format="table",
        data_columns=["CAG"]
    )

    if corncob_df is not None:
        corncob_df.to_hdf(
            store, 
            "/stats/cag/corncob", 
            format="fixed"
        )

"""

}


// Repack an HDF5 file
process repackHDF {

    container "${container__pandas}"
    tag "Compress HDF store"
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