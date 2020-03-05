
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
    take:
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
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

cag_csv = "${cag_csv}"
gene_abund_feather = "${gene_abund_feather}"
cag_abund_feather = "${cag_abund_feather}"
famli_json_list = "${famli_json_list}".split(" ")
readcount_csv = "${readcount_csv}"
manifest_csv = "${manifest_csv}"
suffix=".json.gz"

# Keep a list of descriptive statistics
summary_dict = dict([
    ("formula", "${params.formula}"),
    ("distance_threshold", "${params.distance_threshold}"),
    ("distance_threshold", "${params.distance_threshold}"),
    ("linkage_type", "${params.linkage_type}"),
    ("sd_mean_cutoff", ${params.sd_mean_cutoff}),
    ("min_identity", ${params.min_identity}),
    ("min_coverage", ${params.min_coverage}),
    ("dmnd_min_identity", ${params.dmnd_min_identity}),
    ("dmnd_min_coverage", ${params.dmnd_min_coverage})
])

# Function to read in the FAMLI output
def read_famli_json(fp):

    # Read in the JSON as a DataFrame
    df = pd.DataFrame(
        json.load(
            gzip.open(
                fp, "rt"
            )
        )
    )

    assert df.shape[0] == df.dropna().shape[0], "Missing values for %s" % (sample_name)

    return df

# Perform different types of ordination and include them in the output HDF
def run_ordination(abund_df, ordination_type):

    if ordination_type == "pca":
        return run_pca(abund_df)
    elif ordination_type == "tsne":
        return run_tsne(abund_df)
    else:
        assert "No code available to run ordination type: %s" % ordination_type

# Dimensionality reduction with PCA
def run_pca(abund_df):
    
    # Initialize the PCA object
    pca = PCA()
    
    # Fit to the data
    pca.fit(abund_df)

    # Make an output DataFrame
    return pd.DataFrame(
        pca.transform(
            abund_df
        ),
        index=abund_df.index.values,
        columns=[
            "PC%d (%s%s)" % (
                ix + 1,
                round(100 * r, 1) if r > 0.01 else "%.1E" % (100 * r),
                '%'
            )
            for ix, r in enumerate(
                pca.explained_variance_ratio_
            )
        ]
    ).reset_index(
    ).rename(columns=dict([('index', 'specimen')]))

# Dimensionality reduction with t-SNE
def run_tsne(abund_df, n_components=2):
    
    # Initialize the TSNE object
    tsne = TSNE(
        n_components=n_components
    )

    # Make an output DataFrame with the transformed data
    return pd.DataFrame(
        tsne.fit_transform(
            abund_df
        ),
        index=abund_df.index.values,
        columns=[
            "t-SNE %d" % (
                ix + 1
            )
            for ix in range(n_components)
        ]
    ).reset_index(
    ).rename(columns=dict([('index', 'specimen')]))

# Open a connection to the HDF5
with pd.HDFStore("${params.output_prefix}.full.hdf5", "w") as store:

    # Keep track of the total number of aligned reads for each sample
    aligned_reads_dict = dict()

    # Read in the complete set of FAMLI results
    for fp in famli_json_list:

        # Get the sample name from the file name
        assert fp.endswith(suffix)
        sample_name = fp[:-len(suffix)]

        # Read from JSON for this sample
        print("Reading gene information for %s from %s" % (sample_name, fp))
        df = read_famli_json(fp)
        
        print("Saving to HDF")
        df.to_hdf(
            store, "/abund/gene/long/%s" % sample_name
        )

        # Record the total number of aligned reads for this sample
        aligned_reads_dict[sample_name] = df["nreads"].sum()
    
    # Read in the summary of the number of reads across all samples
    readcount_df = pd.read_csv(readcount_csv)

    # Add some descriptive statistics
    summary_dict["num_samples"] = readcount_df.shape[0]
    summary_dict["total_reads"] = readcount_df["n_reads"].sum()

    # Add the number of aligned reads per sample
    readcount_df["aligned_reads"] = readcount_df["specimen"].apply(
        aligned_reads_dict.get
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
    # Add some descriptive statistics
    summary_dict["num_cags"] = cag_abund_df.shape[0]


    # Write to HDF5
    cag_abund_df.to_hdf(store, "/abund/cag/wide")

    # Perform ordination, both PCA and t-SNE and write to store
    for ordination_type in ["pca", "tsne"]:

        # The ordination methods expect samples to be in rows
        # and so we need to rotate the /abund/cags/wide object
        run_ordination(
            cag_abund_df.set_index("CAG").T,
            ordination_type
        ).to_hdf(
            store,
            "/ordination/%s" % ordination_type
        )

    # Read in the table with the gene-level abundances across all samples
    gene_abund_df = pd.read_feather(
        gene_abund_feather
    )
    print(
        "Read in abundances for %d genes across %d samples" %
        (gene_abund_df.shape[0], gene_abund_df.shape[1] - 1)
    )
    # Add some descriptive statistics
    summary_dict["num_genes"] = gene_abund_df.shape[0]


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

    # Make a summary table describing each CAG with size, mean_abundance, and prevalence
    # Save it in the HDF5 as "/annot/cag/all"
    cag_abund_df.set_index("CAG", inplace=True)
    cag_summary_df = pd.DataFrame(dict([
        ("size", cag_df["CAG"].value_counts()),
        ("prevalence", (cag_abund_df > 0).mean(axis=1)),
        ("mean_abundance", cag_abund_df.mean(axis=1))
    ])).reset_index(
    ).rename(
        columns=dict([("index", "CAG")])
    ).to_hdf(
        store,
        "/annot/cag/all"
    )

    # Write out the descriptive statistics
    pd.DataFrame([
        dict([
            ("variable", k),
            ("value", v)
        ])
        for k, v in summary_dict.items()
    ]).to_hdf(
        store,
        "/summary/experiment"
    )

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
    sep = "\\t",
    compression = "gzip"
)
print(
    "Read in a table matching %d alleles against %d genes" % 
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

# Make a wide version of the table with just the mu. values
corncob_wide = corncob_df.loc[
    corncob_df["parameter"].apply(
        lambda s: s.startswith("mu.")
    )
].pivot_table(
    index = ["CAG", "parameter"],
    columns = "type",
    values = "value"
).reset_index(
).apply(
    lambda v: v.apply(lambda s: s.replace("mu.", "")) if v.name == "parameter" else v
)


# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Write corncob results to HDF5
    corncob_df.to_hdf(store, "/stats/cag/corncob")

    # Write corncob results to HDF5 in wide format
    corncob_wide.to_hdf(store, "/stats/cag/corncob_wide")

"""

}

process addEggnogResults {
    tag "Add functional predictions to HDF"
    container "${container__experiment_collection}"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
        path results_hdf
        path "genes.emapper.annotations.*.gz"

    output:
        path "${results_hdf}"

"""
#!/usr/bin/env python3

import os
import pandas as pd

# Get the list of files with eggNOG results
eggnog_csv_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("genes.emapper.annotations")
]

# Make sure that every eggNOG file ends with '.gz'
# If not, this is a sign that the file was staged incorrectly
for fp in eggnog_csv_list:
    assert fp.endswith('.gz'), "Unexpected: %s" % fp

print("Found %d files with eggNOG results" % len(eggnog_csv_list))

# Read in the eggNOG-mapper results
eggnog_df = pd.concat([
    pd.read_csv(
        fp, 
        header=3, 
        sep="\\t"
    ).rename(columns=dict([("#query_name", "query_name")]))
    for fp in eggnog_csv_list
])

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
        path "genes.tax.aln.*.gz"
        path taxonomy_csv

    output:
        path "${results_hdf}"

"""
#!/usr/bin/env python3

import os
import pandas as pd

diamond_tax_csv_list = [
    fp
    for fp in os.listdir(".")
    if fp.startswith("genes.tax.aln")
]
for fp in diamond_tax_csv_list:
    assert fp.endswith(".gz")

print("Found %d taxonomy results CSVs to import" % len(diamond_tax_csv_list))

# Read in the DIAMOND-tax results
tax_df = pd.concat([
    pd.read_csv(
        fp, 
        sep="\\t", 
        header=None
    ).rename(
        columns=dict([
            (0, "gene"), 
            (1, "tax_id"), 
            (2, "evalue")
        ])
    )
    for fp in diamond_tax_csv_list
])

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
    print("Reading the manifest")
    manifest_df = pd.read_hdf(store, "/manifest")

    # Read in the summary of the entire experiment
    print("Reading the summary")
    exp_summary_df = pd.read_hdf(store, "/summary/experiment")

    # Table with the abundance of every CAG in every sample
    print("Reading the CAG abundance")
    cag_abund_df = pd.read_hdf(store, "/abund/cag/wide")

    # Table with PCA for every sample, based on CAG abundances
    print("Reading the PCA")
    pca_df = pd.read_hdf(store, "/ordination/pca")

    # Table with t-SNE for every sample, based on CAG abundances
    print("Reading the t-SNE")
    tsne_df = pd.read_hdf(store, "/ordination/tsne")

    # Table with gene annotations (including CAG membership)
    print("Reading the gene annotations")
    gene_annot_df = pd.read_hdf(store, "/annot/gene/all")

    # Table summarizing the size, prevalence, and mean abundance of every CAG
    print("Reading the CAG annotations")
    cag_summary_df = pd.read_hdf(store, "/annot/cag/all")

    # Number of read pairs
    print("Reading the readcounts")
    readcount_df = pd.read_hdf(store, "/summary/readcount")

    # Check to see if we have any corncob results
    if "/stats/cag/corncob" in store:

        # Read in the table of corncob results
        print("Reading the statistical analysis results")
        corncob_df = pd.read_hdf(store, "/stats/cag/corncob")

        # Read in the table of corncob results in wide format
        print("Reading the statistical analysis results in wide format")
        corncob_wide = pd.read_hdf(store, "/stats/cag/corncob_wide")

    else:
        corncob_df = None
        corncob_wide = None

    # Check to see if the taxonomy was included in the outputs
    if "/ref/taxonomy" in store:

        # Read in the taxonomy
        taxonomy_df = pd.read_hdf(store, "/ref/taxonomy")

    else:
        taxonomy_df = None

# Now write to the new HDF store
with pd.HDFStore("${params.output_prefix}.summary.hdf5", "w") as store:

    # Write out each table, while indexing the appropriate columns
    print("Writing the manifest")
    manifest_df.to_hdf(
        store,
        "/manifest"
    )
    print("Writing the summary")
    exp_summary_df.to_hdf(
        store,
        "/summary/experiment"
    )
    print("Writing the readcounts")
    readcount_df.to_hdf(
        store,
        "/summary/readcount"
    )
    print("Writing the CAG abundance")
    cag_abund_df.to_hdf(
        store, 
        "/abund/cag/wide", 
        format="table",
        data_columns=["CAG"]
    )
    print("Writing the PCA")
    pca_df.to_hdf(
        store,
        "/ordination/pca",
        format="fixed"
    )
    print("Writing the t-SNE")
    tsne_df.to_hdf(
        store,
        "/ordination/tsne",
        format="fixed"
    )
    print("Writing the gene annotations")
    gene_annot_df.to_hdf(
        store, 
        "/annot/gene/all", 
        format="table",
        data_columns=["CAG"]
    )
    print("Writing the CAG annotations")
    cag_summary_df.to_hdf(
        store,
        "/annot/cag/all",
        format="fixed"
    )

    if corncob_df is not None:
        print("Writing the statistical analysis results")
        corncob_df.to_hdf(
            store, 
            "/stats/cag/corncob", 
            format="fixed"
        )
        print("Writing the statistical analysis results in wide format")
        corncob_wide.to_hdf(
            store, 
            "/stats/cag/corncob_wide", 
            format="fixed"
        )

    if taxonomy_df is not None:
        print("Writing the taxonomy table")

        taxonomy_df.to_hdf(
            store,
            "/ref/taxonomy",
            format = "fixed"
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