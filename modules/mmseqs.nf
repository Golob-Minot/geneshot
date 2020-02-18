container__pandas = "quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed"

process mmseqs {
    tag "Cluster genes with similar sequences"
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'mem_medium'
    errorStrategy 'retry'
    
    input:
    file "input.genes.fasta.gz"
    
    output:
    file "output.genes.fasta.gz"
    file "output.genes.alleles.tsv.gz"
    
"""
#!/bin/bash

set -e

# Make the MMSeqs2 database
mmseqs createdb input.genes.fasta.gz db

# Cluster the protein sequences
mmseqs linclust db cluster_db ./ \
    --min-seq-id ${params.min_identity / 100} \
    --max-seqs 100000 \
    -c ${params.min_coverage / 100}

# Make TSV output for clustering
mmseqs createtsv db db cluster_db output.genes.alleles.tsv
echo "Compressing"
gzip output.genes.alleles.tsv

# Get the representative sequences
mmseqs result2repseq db cluster_db genes
mmseqs result2flat db db genes output.genes.fasta --use-fasta-header
echo "Compressing"
gzip output.genes.fasta

echo "Done"
"""
}

process joinGeneClusters {
    tag "Combine multiple rounds of gene clustering"
    container "${container__pandas}"
    label 'mem_medium'
    // errorStrategy 'retry'

    input:
    file "genes.alleles.*.tsv.gz"

    output:
    file "genes.alleles.tsv.gz"

"""
#!/usr/bin/env python3

import pandas as pd
import os

# Input files will be staged as genes.alleles.1.tsv.gz, genes.alleles.2.tsv.gz, etc.

# Read in all of the tables
cluster_list = []

# Iterate over all of the files in the working directory
for fp in os.listdir("."):

    # Make sure the file has the right name
    if fp.startswith("genes.alleles.") is False:
        continue
    assert fp.endswith(".tsv.gz"), "%s is unexpected" % (fp)

    # Read in the table
    print("Reading in %s" % fp)
    cluster_list.append(
        pd.read_csv(
            fp, 
            sep="\\t", 
            compression="gzip", 
            header=None,
            names=['gene', 'allele'], 
            usecols=[0, 1]
        )
    )
    print(cluster_list[-1].head())

# Sort the clusters by size, descending
cluster_list.sort(
    key=lambda df: df.shape[0],
    reverse=True
)

# Logging and sanity checking
for ix, df in enumerate(cluster_list):
    print("Cluster %d has %d alleles grouped into %d genes" % (ix + 1, df.shape[0], df["gene"].unique().shape[0]))

    # The number of alleles in this cluster is the number of genes in the previous cluster
    if ix > 0:
        assert df.shape[0] == cluster_list[ix - 1]["gene"].unique().shape[0]

# Make a master set of clusters, starting with the first
cluster_table = cluster_list[0]
print("Starting with a master table of %d alleles grouped into %d genes" % (cluster_table.shape[0], cluster_table["gene"].unique().shape[0]))

# Now iteratively add in the subsequent clusters
for ix, df in enumerate(cluster_list):
    if ix == 0:
        continue

    # Replace the 'genes' in the master list with the 'genes' found in subsequent rounds
    cluster_table.replace(
        to_replace = dict([
            ('gene', df.query("gene != allele").set_index("allele")["gene"].to_dict())
        ]),
        inplace=True
    )

    print("Master table now has %d alleles grouped into %d genes" % (cluster_table.shape[0], cluster_table["gene"].unique().shape[0]))

# Write out the final set of clusters
print("Writing out to genes.alleles.tsv.gz")
cluster_table.to_csv("genes.alleles.tsv.gz", index=None, header=False, sep="\\t")
print("Done")
"""
}