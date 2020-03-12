
// Container versions
container__fastatools = "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.2"
container__ubuntu = "ubuntu:18.04"
container__experiment_collection = "quay.io/fhcrc-microbiome/experiment-collection@sha256:fae756a380a3d3335241b68251942a8ed0bf1ae31a33a882a430085b492e44fe"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed"

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
        path famli_json_list
        path readcount_csv
        path manifest_csv
        path specimen_gene_count_csv
        path breakaway_csv

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
specimen_gene_count_csv = "${specimen_gene_count_csv}"
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

    # Keep track of the length of each gene
    gene_length_dict = {}

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

        # Add the gene lengths
        for _, r in df.iterrows():
            gene_length_dict[r["id"]] = r["length"]

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
    # This is being called 'cag_df', but it's really a table of CAG annotations per-gene,
    # so there is one row per gene.
    cag_df = pd.read_csv(cag_csv)

    print(
        "Read in CAG assignments for %d genes across %d CAGs" % 
        (cag_df.shape[0], cag_df["CAG"].unique().shape[0])
    )

    # Write to HDF5
    cag_df.to_hdf(store, "/annot/gene/cag")

    # Calculate prevalence and abundance information for each gene
    gene_abund_df.set_index("index", inplace=True)
    gene_abundance = gene_abund_df.mean(axis=1)
    gene_prevalence = (gene_abund_df > 0).mean(axis=1)

    # Add that information on the gene abundance and prevalence to this gene summary table
    cag_df = cag_df.assign(
        prevalence = cag_df["gene"].apply(gene_prevalence.get),
        abundance = cag_df["gene"].apply(gene_abundance.get),
        length = cag_df["gene"].apply(gene_length_dict.get)
    )

    # Now create the `/annot/gene/all` table, which will be added to later
    # with the taxonomic and functional annotations, if those are performed
    cag_df.to_hdf(store, "/annot/gene/all")

    # Make a summary table describing each CAG with size, mean_abundance, and prevalence
    # Save it in the HDF5 as "/annot/cag/all"
    cag_abund_df.set_index("CAG", inplace=True)
    cag_summary_df = pd.DataFrame(dict([
        ("size", cag_df["CAG"].value_counts()),
        ("prevalence", (cag_abund_df > 0).mean(axis=1)),
        ("mean_abundance", cag_abund_df.mean(axis=1)),
        ("std_abundance", cag_abund_df.std(axis=1))
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

# Count up the number of genes assembled per-specimen
n_genes_assembled_per_specimen = pd.DataFrame(
    dict(
        [
            (
                "n_genes_assembled",
                allele_assembly["specimen"].value_counts()
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

    # Write the summary of the number of genes assembled per sample
    n_genes_assembled_per_specimen.to_hdf(store, "/summary/genes_assembled")

    print(pd.read_hdf(
        store,
        "/summary/all"
    ).head())

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

# Keep track of data objects in a dict
# This should include the key, the format to write the HDF, and any data_columns if needed
# This is the list of things which will be written to the summary HDF5,
# noting those objects which are optional or required
data_objects = [
    dict([
        ("key", "/manifest"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "manifest describing experiment design")
    ]),
    dict([
        ("key", "/summary/experiment"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "informatics parameter summary")
    ]),
    dict([
        ("key", "/summary/readcount"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "number of read pairs")
    ]),
    dict([
        ("key", "/summary/genes_assembled"),
        ("optional", True),
        ("format", "fixed"),
        ("comment", "number of genes identified by de novo assembly per sample")
    ]),
    dict([
        ("key", "/summary/genes_aligned"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "number of genes identified by alignment per sample")
    ]),
    dict([
        ("key", "/summary/breakaway"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "gene richness estimated by breakaway")
    ]),
    dict([
        ("key", "/summary/all"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "combined summary metrics for all specimens")
    ]),
    dict([
        ("key", "/abund/cag/wide"),
        ("optional", False),
        ("format", "table"),
        ("data_columns", ["CAG"]),
        ("comment", "abundance of every CAG in every sample")
    ]),
    dict([
        ("key", "/ordination/pca"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "PCA for every sample, based on CAG abundances")
    ]),
    dict([
        ("key", "/ordination/tsne"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "t-SNE for every sample, based on CAG abundances")
    ]),
    dict([
        ("key", "/annot/gene/all"),
        ("optional", False),
        ("format", "table"),
        ("data_columns", ["CAG"]),
        ("comment", "table with gene annotations")
    ]),
    dict([
        ("key", "/annot/cag/all"),
        ("optional", False),
        ("format", "fixed"),
        ("comment", "table with CAG annotations")
    ]),
    dict([
        ("key", "stats/cag/corncob"),
        ("optional", True),
        ("format", "fixed"),
        ("comment", "table of corncob results")
    ]),
    dict([
        ("key", "/stats/cag/corncob_wide"),
        ("optional", True),
        ("format", "fixed"),
        ("comment", "table of corncob results in wide format")
    ]),
    dict([
        ("key", "/ref/taxonomy"),
        ("optional", True),
        ("format", "fixed"),
        ("comment", "taxonomy table")
    ]),
]

# Open a connection to the input store
input_store = pd.HDFStore("${full_results_hdf}", "r")

# Open a connection to the output store
output_store = pd.HDFStore("${params.output_prefix}.summary.hdf5", "w")

# Iterate over each of the data objects, copying from the input to the output
for d in data_objects:
    # Make sure the key is present
    if d["key"] in input_store:
        # Copy the table to the output store
        print("Reading in %s and writing to the summary HDF5" % d["comment"])
        if d["format"] == "fixed":
            pd.read_hdf(
                input_store, 
                d["key"]
            ).to_hdf(
                output_store,
                d["key"],
                format = d["format"]
            )
        else:
            assert d["format"] == "table", d
            assert "data_columns" in d
            pd.read_hdf(
                input_store, 
                d["key"]
            ).to_hdf(
                output_store,
                d["key"],
                format = d["format"],
                data_columns = d["data_columns"]
            )
    else:
        assert d["optional"], "Could not find expected key %s containing %s" % (d["key"], d["comment"])
        print("Did not find the key %s in the HDF5" % d["key"])

# Close the stores
input_store.close()
output_store.close()
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