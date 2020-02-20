container__pandas = "quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed"

process linclust {
    tag "Cluster genes with similar sequences"
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'mem_medium'
    errorStrategy 'retry'
    
    input:
    file "input.genes.*.fasta.gz"
    
    output:
    file "output.genes.fasta.gz"
    
"""
#!/bin/bash

set -e

# Combine input files
echo "Combining input files"
cat input.genes.* > input.genes.fasta.gz

# Make the MMSeqs2 database
echo "Running linclust"
mmseqs createdb input.genes.fasta.gz db

# Cluster the protein sequences
mmseqs linclust db cluster_db ./ \
    --min-seq-id ${params.min_identity / 100} \
    --max-seqs 100000 \
    -c ${params.min_coverage / 100}

# Get the representative sequences
mmseqs result2repseq db cluster_db genes
mmseqs result2flat db db genes output.genes.fasta --use-fasta-header
echo "Compressing"
gzip output.genes.fasta

echo "Done"
"""
}
