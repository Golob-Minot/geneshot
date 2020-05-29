
// Container versions
container__fastatools = "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.2"
container__ubuntu = "ubuntu:18.04"
container__experiment_collection = "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

// Default parameters
params.fdr_method = "fdr_bh"

// Function to read in a CSV and return a Channel
def read_manifest(manifest_file){
    manifest_file.splitCsv(
        header: true, 
        sep: ","
    ).branch{
        valid_paired_indexed:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty()) && (it.R2 != null ) && (it.R2 != "") && (!file(it.R2).isEmpty()) && (it.I1 != null ) && (it.I1 != "" ) && (!file(it.I1).isEmpty()) && (it.I2 != null ) && (it.I2 != "") && (!file(it.I2).isEmpty())
        valid_paired:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty()) && (it.R2 != null ) && (it.R2 != "") && (!file(it.R2).isEmpty())
        valid_unpaired:  (it.specimen != null) && (it.R1 != null ) && (it.R1 != "" ) && (!file(it.R1).isEmpty())
        other: true
    }
}

workflow combineReads {
    take:

        fastq_ch

    main:

        fastq_ch.branch {  // Split up the samples which have multiple FASTQ files
            single: it[1].size() == 1
            multiple: it[1].size() > 1
        }.set {
            grouped_fastq
        }

        joinFASTQ(
            grouped_fastq.multiple
        )

    emit:
        grouped_fastq.single.map {
            r -> [r[0], r[1][0], r[2][0]]
        }.mix(
            joinFASTQ.out
        )

}

process joinFASTQ {
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

(( \$(gunzip -c "${sample}.R1.fastq.gz" | head | wc -l) > 1 ))
(( \$(gunzip -c "${sample}.R2.fastq.gz" | head | wc -l) > 1 ))

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
        path readcount_csv
        path manifest_csv
        path specimen_gene_count_csv
        path specimen_reads_aligned_csv
        path gene_length_csv
        path breakaway_csv

    output:
        path "${params.output_prefix}.results.hdf5"

"""
#!/usr/bin/env python3

import gzip
import json
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from scipy.spatial.distance import pdist, squareform
from scipy.stats.mstats import gmean
from scipy.stats import entropy

cag_csv = "${cag_csv}"
gene_abund_feather = "${gene_abund_feather}"
cag_abund_feather = "${cag_abund_feather}"
readcount_csv = "${readcount_csv}"
specimen_gene_count_csv = "${specimen_gene_count_csv}"
specimen_reads_aligned_csv = "${specimen_reads_aligned_csv}"
gene_length_csv = "${gene_length_csv}"
breakaway_csv = "${breakaway_csv}"
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

# Calculate pairwise distances (input has specimens in rows)
def calc_pdist(abund_df, metric="euclidean"):
    print("Calculating pairwise %s distances over %d samples and %d CAGs" % (metric, abund_df.shape[0], abund_df.shape[1]))
    
    # Make an output DataFrame with the pairwise distances
    return pd.DataFrame(
        squareform(
            pdist(
                abund_df if metric != "aitchison" else clr_transform(abund_df),
                metric=metric if metric != "aitchison" else "euclidean"
            )
        ),
        index=abund_df.index.values,
        columns=abund_df.index.values
    ).reset_index(
    ).rename(columns=dict([('index', 'specimen')]))

# Calculate CLR-transformed abundances (input has specimens in rows)
def clr_transform(abund_df):
    print("Calculating CLR over %d samples and %d CAGs" % (abund_df.shape[0], abund_df.shape[1]))

    # Make the input
    specimen_col_abund_df = abund_df.T

    # Calculate the geometric mean for every specimen
    specimen_gmean = specimen_col_abund_df.apply(
        lambda c: gmean(c.loc[c > 0])
    )

    # Find the global minimum value
    min_abund = specimen_col_abund_df.apply(
        lambda c: c.loc[c > 0].min()
    ).min()

    print("Minimum value is %d" % min_abund)

    # Fill in the minimum value
    specimen_col_abund_df = specimen_col_abund_df.clip(
        lower=min_abund
    )

    # Divide by the geometric mean
    specimen_col_abund_df = specimen_col_abund_df / specimen_gmean
    
    # Take the log and return
    return specimen_col_abund_df.applymap(np.log10).T

# Read in the table with gene lengths
gene_length_df = pd.read_csv(gene_length_csv).set_index("gene")["length"]

# Open a connection to the HDF5
with pd.HDFStore("${params.output_prefix}.results.hdf5", "w") as store:

    # Read in the summary of the number of reads across all samples
    readcount_df = pd.read_csv(readcount_csv)

    # Read in the number of aligned reads per sample
    aligned_reads_dict = pd.read_csv(
        specimen_reads_aligned_csv
    ).set_index(
        "specimen"
    )["n_reads_aligned"]

    # Add some descriptive statistics
    summary_dict["num_samples"] = readcount_df.shape[0]
    summary_dict["total_reads"] = readcount_df["n_reads"].sum()
    summary_dict["aligned_reads"] = aligned_reads_dict.sum()

    # Add the number of aligned reads per sample
    readcount_df["aligned_reads"] = readcount_df["specimen"].apply(
        aligned_reads_dict.get
    )

    # Write to HDF5
    readcount_df.to_hdf(store, "/summary/readcount")

    # Read in the summary of the number of genes detected by alignment per sample
    specimen_gene_count_df = pd.read_csv(specimen_gene_count_csv)

    # Write to HDF5
    specimen_gene_count_df.to_hdf(store, "/summary/genes_aligned")

    # Read in the gene richness estimate generated with breakaway
    breakaway_df = pd.read_csv(breakaway_csv)

    # Write to HDF5
    breakaway_df.to_hdf(store, "/summary/breakaway")

    # Write out a combined table
    pd.concat(
        [
            readcount_df.set_index("specimen"),
            specimen_gene_count_df.set_index("specimen"),
            breakaway_df.set_index("specimen")
        ], 
        axis = 1, 
        sort = True
    ).reset_index(
    ).rename(
        columns = dict([
            ("index", "specimen")
        ])
    ).to_hdf(
        store,
        "/summary/all"
    )

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
    cag_abund_df.to_hdf(
        store, 
        "/abund/cag/wide",
        format = "table",
        data_columns = ["CAG"],
        complevel = 5
    )

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

    # Calculate pairwise distances and write to store
    for metric in ["euclidean", "braycurtis", "jaccard", "aitchison"]:

        # Calculate pairwise distances and write to the store
        calc_pdist(
            cag_abund_df.set_index("CAG").T,
            metric
        ).to_hdf(
            store,
            "/distances/%s" % metric
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
    # This is being called 'cag_df', but it's really a table of CAG annotations per-gene,
    # so there is one row per gene.
    cag_df = pd.read_csv(cag_csv)

    print(
        "Read in CAG assignments for %d genes across %d CAGs" % 
        (cag_df.shape[0], cag_df["CAG"].unique().shape[0])
    )

    # Write to HDF5
    cag_df.to_hdf(
        store, 
        "/annot/gene/cag",
        format = "table",
        data_columns = ["gene", "CAG"]
    )

    # Calculate prevalence and abundance information for each gene
    gene_abund_df.set_index("index", inplace=True)
    gene_abundance = gene_abund_df.mean(axis=1)
    gene_prevalence = (gene_abund_df > 0).mean(axis=1)

    # Add that information on the gene abundance and prevalence to this gene summary table
    cag_df = cag_df.assign(
        prevalence = cag_df["gene"].apply(gene_prevalence.get),
        abundance = cag_df["gene"].apply(gene_abundance.get),
        length = cag_df["gene"].apply(gene_length_df.get)
    )

    # Now create the `/annot/gene/all` table, which will be added to later
    # with the taxonomic and functional annotations, if those are performed
    cag_df.to_hdf(
        store, 
        "/annot/gene/all",
        format = "table",
        data_columns = ["gene", "CAG"]
    )

    # Make a summary table describing each CAG with size, mean_abundance, and prevalence
    # Save it in the HDF5 as "/annot/cag/all"
    cag_abund_df.set_index("CAG", inplace=True)
    cag_summary_df = pd.DataFrame(dict([
        ("size", cag_df["CAG"].value_counts()),
        ("prevalence", (cag_abund_df > 0).mean(axis=1)),
        ("mean_abundance", cag_abund_df.mean(axis=1)),
        ("std_abundance", cag_abund_df.std(axis=1)),
        ("entropy", cag_abund_df.apply(entropy, axis=1)),
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

    # Read in the manifest provided by the user
    manifest_df = pd.read_csv(
        "${manifest_csv}"
    )
    # Drop the columns with paths to reads
    for col_name in ["R1", "R2", "I1", "I2"]:
        if col_name in manifest_df.columns.values:
            manifest_df = manifest_df.drop(columns=col_name)

    # Drop duplicated rows (multiple read pairs from the same set of samples)
    manifest_df = manifest_df.drop_duplicates()

    # Write out the manifest provided by the user
    manifest_df.to_hdf(store, "/manifest")

"""

}

process addGeneAssembly{
    tag "Add gene assembly data to HDF"
    container "${container__experiment_collection}"
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
        path results_hdf
        path detailed_hdf
        path allele_assembly_csv_list

    output:
        path "${results_hdf}"
        path "${detailed_hdf}"

"""
#!/usr/bin/env python3

import pandas as pd
import os

# Count up the number of genes assembled per-specimen
n_genes_assembled_per_specimen = dict()

# Open a connection to the detailed HDF5 output store
detailed_store = pd.HDFStore("${detailed_hdf}", "a")

# Read in the summary of allele assembly for each sample
for fp in "${allele_assembly_csv_list}".split(" "):

    # Make sure that the file path has the expected pattern
    assert fp.endswith(".csv.gz"), fp
    sample_name = fp.replace(".csv.gz", "")
    print("Reading in assembly information for %s" % sample_name)
    assembly_df = pd.read_csv(fp)
    print("Read in %d assembled genes" % assembly_df.shape[0])

    # Save the information on how many genes were assembled in each sample
    n_genes_assembled_per_specimen[sample_name] = assembly_df.shape[0]

    # Write out to the detailed HDF
    assembly_df.to_hdf(
        detailed_store,
        "/abund/allele/assembly/%s" % sample_name
    )

detailed_store.close()

n_genes_assembled_per_specimen = pd.DataFrame(
    dict(
        [
            (
                "n_genes_assembled",
                pd.Series(n_genes_assembled_per_specimen)
            )
        ]
    )
).reset_index(
).rename(
    columns = dict(
        [
            ("index", "specimen")
        ]
    )
)

# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Write the summary of the number of genes assembled per sample
    n_genes_assembled_per_specimen.to_hdf(store, "/summary/genes_assembled")

    # Add the number of genes assembled to the combined summary table
    pd.concat([
        pd.read_hdf(
            store,
            "/summary/all"
        ).set_index(
            "specimen"
        ),
        n_genes_assembled_per_specimen.set_index(
            "specimen"
        )
    ], axis = 1, sort = True).reset_index(
    ).rename(
        columns = dict([
            ("index", "specimen")
        ])
    ).to_hdf(
        store,
        "/summary/all"
    )

"""

}

process addCorncobResults{
    tag "Add statistical analysis to HDF"
    container "${container__pandas}"
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
from statsmodels.stats.multitest import multipletests

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

# Adding the q-value is conditional on p-values being present
if "p_value" in corncob_wide.columns.values:

    # Add the q-value (FDR-BH)
    corncob_wide = corncob_wide.assign(
        q_value = multipletests(corncob_wide.p_value.fillna(1), 0.2, "${params.fdr_method}")[1]
    )

# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Write corncob results to HDF5
    corncob_wide.to_hdf(store, "/stats/cag/corncob")

"""

}

process addMetaPhlAn2Results{
    tag "Add composition analysis to HDF"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'retry'

    input:
        path results_hdf
        path metaphlan_tsv_list

    output:
        path "${results_hdf}"

"""
#!/usr/bin/env python3

import pandas as pd

# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Iterate over each input file
    for fp in "${metaphlan_tsv_list}".split(" "):

        # Make sure that the file has the expected suffix
        assert fp.endswith(".metaphlan2.tsv"), fp

        # Get the specimen name from the file name
        specimen_name = fp.replace(".metaphlan2.tsv", "")

        # Read in from the flat file
        metaphlan_df = pd.read_csv(
            fp,
            skiprows=1,
            sep="\\t"
        )

        print("Read in %d taxa from %s" % (metaphlan_df.shape[0], specimen_name))

        # Write metaphlan results to HDF5
        key = "/composition/metaphlan/%s" % specimen_name
        print("Writing to %s" % key)
        metaphlan_df.to_hdf(store, key)

    print("Closing store")

print("Done")

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
    )
    for fp in eggnog_csv_list
]).rename(
    columns = dict([("#query_name", "query_name")])
).reset_index(
    drop=True
)

def head_and_tail(df):
    print("HEAD:")
    print(df.head())
    print("")
    print("TAIL:")
    print(df.tail())
    print("")

head_and_tail(eggnog_df)

assert 'query_name' in eggnog_df.columns.values
print("Read in %d lines of annotations" % eggnog_df.shape[0])

# Remove those rows with no "seed_eggNOG_ortholog" (which are internal runtime metrics)
eggnog_df = eggnog_df.reindex(
    index=eggnog_df["seed_eggNOG_ortholog"].dropna().index
)
print("Read in annotations for %d genes" % eggnog_df.shape[0])

head_and_tail(eggnog_df)

# Set the index to be the gene name
eggnog_df.set_index("query_name", inplace=True)
print("Set the index on 'query_name'")

head_and_tail(eggnog_df)

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
    eggNOG_ortholog = eggnog_df["seed_eggNOG_ortholog"]
).assign(
    eggNOG_tax = eggnog_df["best_tax_level"]
).assign(
    eggNOG_desc = eggnog_df["eggNOG free text desc."]
).sort_values(
    by="CAG"
).reset_index()

# Open a connection to the HDF5
with pd.HDFStore("${results_hdf}", "a") as store:

    # Write eggnog results to HDF5
    eggnog_df.reset_index().reindex(
        columns = [
            'query_name', 
            'seed_eggNOG_ortholog', 
            'best_tax_level', 
            'GOs', 
            'EC', 
            'KEGG_ko', 
            'best eggNOG OG', 
            'COG Functional cat.', 
            'eggNOG free text desc.'
        ]
    ).to_hdf(
        store,
        "/annot/gene/eggnog",
        format = "fixed",
        dtype=str,
    )
    
    # Write summary annotation table to HDF5
    gene_annot.to_hdf(
        store, 
        "/annot/gene/all",
        format = "table",
        data_columns = ["CAG"],
        dtype=str,
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

print("Opening taxdump tar file")
tar = tarfile.open(
    "${taxdump_tar_gz}", 
    "r:gz", 
    encoding='utf-8'
)

# Read in the table of merged taxids
print("Reading in merged taxids")
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
print("Reading in names")
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
print("Reading in nodes")
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
print("Joining names and nodes")
tax_df = pd.concat([
    names_df.set_index("tax_id"),
    nodes_df.set_index("tax_id")
], axis=1, sort=True)

# Add in the merged taxids
print("Adding in merged tax IDs")
tax_df = pd.concat([
    tax_df,
    tax_df.reindex(
        index=merged_df["old"].values
    ).dropna(
    ).apply(
        lambda v: v.apply(merged_df.set_index("old")["new"].get) if v.name == "tax_id" else v
    )
]).reset_index()
assert tax_df["tax_id"].apply(lambda n: "\\n" not in str(n)).all()

# Write to CSV
print("Writing out final CSV")
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
        "/annot/gene/all",
        format = "table",
        data_columns = ["gene", "CAG"]
    )

    # Write the taxonomy table
    taxonomy_df.to_hdf(
        store,
        "/ref/taxonomy"
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