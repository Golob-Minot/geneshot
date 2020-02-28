#!/usr/bin/env nextflow

// Function which prints help message text
def helpMessage() {
    log.info"""

    Geneshot: A pipeline to robustly identify which alleles (n.e.e peptide coding sequences)
    are present in a microbial community.

    One feature of Geneshot is the annotation of coding sequences on the basis of their
    inferred taxonomic identity (taxid) and biophysical function. This workflow can be used
    to generate those reference databases for use in the primary geneshot workflow.

    Taxonomic identity is inferred using the DIAMOND aligner and the NCBI RefSeq database.
    File name: <PREFIX>.refseq.tax.dmnd

    Biophysical function is inferred using eggNOG.
    File names: <PREFIX>.eggnog.db, <PREFIX>.eggnog_proteins.dmnd

    Usage:

    nextflow run jgolob/geneshot/make_reference_databases.nf --output_folder <>
    
    Required Arguments:
      --output_folder       Folder to place database outputs (default ./db/)
      --output_prefix       Prefix to prepend to database outputs (default DB)

    """.stripIndent()
}

params.output_folder = "./ref/"
params.output_prefix = "DB"


process get_taxonomic_db {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'
    publishDir params.output_folder

    output:
    path "${params.output_prefix}.refseq.tax.dmnd"

"""
set -e

echo "Fetching ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip

echo "Unzipping taxdmp.zip"
unzip taxdmp.zip

echo "Making sure that nodes.dmp was found in taxdmp.zip"
ls -lahtr
[[ -s nodes.dmp ]]

echo "Fetching ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

echo "Fetching the NCBI RefSeq non-redundant proteins"
wget -r -l1 -np "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/" -P ./ -A "bacteria.nonredundant_protein.*.protein.faa.gz"

echo "Combining downloaded files"
for fp in ./ftp.ncbi.nlm.nih.gov/refseq/release/bacteria/bacteria.nonredundant_protein.*.protein.faa.gz; do
    
    cat \$fp >> ref.faa.gz
    
    rm \$fp

done

echo "Indexing with DIAMOND"
diamond \
    makedb \
    --in ref.faa.gz \
    --db ${params.output_prefix}.refseq.tax.dmnd \
    --threads ${task.cpus} \
    --taxonmap prot.accession2taxid.gz \
    --taxonnodes nodes.dmp

"""
}

process get_eggnog_db {
    container "quay.io/biocontainers/eggnog-mapper:2.0.1--py_1"
    label 'io_limited'
    publishDir params.output_folder

    output:
    path "${params.output_prefix}.eggnog.db"
    path "${params.output_prefix}.eggnog_proteins.dmnd"

"""
set -e

download_eggnog_data.py -y --data_dir ./

mv eggnog.db ${params.output_prefix}.eggnog.db
mv eggnog_proteins.dmnd ${params.output_prefix}.eggnog_proteins.dmnd


"""
}