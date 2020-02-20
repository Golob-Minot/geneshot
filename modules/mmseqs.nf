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
    label 'mem_veryhigh'
    errorStrategy 'retry'

    input:
    file "old.alleles.tsv.gz"
    file "new.alleles.tsv.gz"

    output:
    file "merged.alleles.tsv.gz"

"""
#!/usr/bin/env python3

import pandas as pd
import os

# Function to read in table clustering alleles into genes
def read_clusters(fp):
    return pd.read_csv(
        fp, 
        sep="\\t", 
        compression="gzip", 
        header=None,
        names=['gene', 'allele'], 
        usecols=[0, 1]
    ).applymap(
        str
    )

# Read in the tables
old_clusters = read_clusters("old.alleles.tsv.gz")
new_clusters = read_clusters("new.alleles.tsv.gz")

# Logging and sanity checking
msg = "The %s set of clusters have %d alleles grouped into %d genes"
print(msg % ("older", old_clusters.shape[0], old_clusters["gene"].unique().shape[0]))
print(msg % ("newer", new_clusters.shape[0], new_clusters["gene"].unique().shape[0]))

# The number of alleles in the newer cluster is the number of genes in the previous cluster
assert new_clusters.shape[0] == old_clusters["gene"].unique().shape[0]

# Replace the 'genes' in the older clusters with the 'genes' found in the newer clusters
if new_clusters.query("gene != allele").shape[0] > 0:
    old_clusters.replace(
        to_replace = dict([
            ('gene', new_clusters.query("gene != allele").set_index("allele")["gene"].to_dict())
        ]),
        inplace=True
    )
else:
    print("Newer clusters did not add any new information, skipping")

print("Combined table now has %d alleles grouped into %d genes" % (old_clusters.shape[0], old_clusters["gene"].unique().shape[0]))

# Write out the final set of clusters
print("Writing out to merged.alleles.tsv.gz")
old_clusters.to_csv("merged.alleles.tsv.gz", index=None, header=False, sep="\\t")
print("Done")
"""
}