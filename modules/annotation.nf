#!/usr/bin/env nextflow

// Processes to perform de novo assembly and annotate those assembled sequences
nextflow.enable.dsl=2


// Using DSL-2
nextflow.enable.dsl=2

// Parameters
params.gene_fasta = false
params.output_hdf = false
params.output_folder = false
params.help = false
params.min_coverage = 50
params.min_identity = 90
params.taxonomic_dmnd = false
params.ncbi_taxdump = "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
params.eggnog_db = false
params.eggnog_dmnd = false
params.noannot = false


process shard_genes {
    tag "Split the gene catalog into smaller shards"
    container "ubuntu:18.04"
    label 'mem_medium'
    errorStrategy 'finish'
    
    input:
    file fasta_gz
    
    output:
    file "genes.shard.*.fasta.gz"
    
"""
#!/bin/bash

set -e

split --additional-suffix .fasta -l 1000000 <(gunzip -c ${fasta_gz}) genes.shard.

gzip genes.shard.*.fasta
"""
}

// Break out the annotation into a separate workflow
workflow Annotation_wf {
    take:
    gene_fasta

    main:

    // Split up the gene catalog into shards for more efficient processing
    shard_genes(
        gene_fasta
    )

    // Determine whether or not to run the eggNOG annotation based
    // on --noannot and --eggnog_db / --eggnog_dmnd
    run_eggnog = false
    if ( params.noannot == false ) {
        if ( params.eggnog_db && params.eggnog_dmnd ) {
            if ( !file(params.eggnog_db).isEmpty() && !file(params.eggnog_dmnd).isEmpty() ){
                run_eggnog = true
            }
        }
    }

    // Annotate the clustered genes with eggNOG
    if ( run_eggnog ){
        eggnog(
            shard_genes.out.flatten(),
            file(params.eggnog_db),
            file(params.eggnog_dmnd)
        )
        eggnog_tsv = eggnog.out.collect()
    } else {
        eggnog_tsv = false
    }

    // Determine whether or not to run the taxnomic annotation based
    // on --noannot and --taxonomic_dmnd
    run_tax = false
    if ( params.noannot == false ) {
        if ( params.taxonomic_dmnd ) {
            if ( !file(params.taxonomic_dmnd).isEmpty() ){
                run_tax = true
            }
        }
    }

    // Annotate the clustered genes with DIAMOND for taxonomic ID
    if ( run_tax ) {
        diamond_tax(
            shard_genes.out.flatten(),
            file(params.taxonomic_dmnd)
        )
        tax_tsv = diamond_tax.out.collect()
        join_tax(
            tax_tsv
        )
    } else {
        tax_tsv = false
    }

    emit:

        tax_tsv = tax_tsv
        eggnog_tsv = eggnog_tsv

}



process diamond_tax {
    tag "Annotate genes by taxonomy"
    container "quay.io/fhcrc-microbiome/famli:v1.5"
    label 'mem_veryhigh'

    input:
    file query
    file diamond_tax_db
    
    output:
    file "genes.tax.aln.gz"

    
"""
set -e

diamond \
    blastp \
    --db ${diamond_tax_db} \
    --query ${query} \
    --out genes.tax.aln.gz \
    --outfmt 102 \
    --id ${params.min_identity} \
    --top ${100 - params.min_identity} \
    --block-size ${task.memory.toMega() / (1024 * 6)} \
    --threads ${task.cpus} \
    --compress 1

rm ${diamond_tax_db}
"""

}

process join_tax {
    tag "Concatenate taxonomy annotation files"
    container "ubuntu:18.04"
    label 'mem_medium'
    publishDir "${params.output_folder}/annot/", mode: "copy"

    input:
    file "genes.tax.aln.*.gz"
    
    output:
    file "genes.tax.aln.gz"

    
"""
set -e

for fp in genes.tax.aln.*.gz; do

    cat \$fp
    rm \$fp

done > genes.tax.aln.gz
"""

}


process eggnog {
    tag "Annotate genes by predicted function"
    container "quay.io/biocontainers/eggnog-mapper:2.1.12--pyhdfd78af_0"
    label 'mem_veryhigh'
    
    input:
    path query
    path eggnog_db
    path eggnog_dmnd

    output:
    path "genes.emapper.annotations.gz"

    
    """
set -e

mkdir data
mkdir TEMP
mkdir SCRATCH

mv ${eggnog_db} data/eggnog.db
mv ${eggnog_dmnd} data/eggnog_proteins.dmnd

emapper.py \
    -i ${query} \
    --itype proteins \
    --output genes \
    -m "diamond" \
    --cpu ${task.cpus} \
    --data_dir data/ \
    --scratch_dir SCRATCH/ \
    --temp_dir TEMP/

gzip genes.emapper.annotations
    
    """

}

// Function which prints help message text
def helpMessage() {
    log.info"""
    Add functional and/or taxonomic annotations to a geneshot output file.

    Usage:

    nextflow run Golob-Minot/geneshot/modules/annotation <ARGUMENTS>
    
    Options:
      --gene_fasta          Location for input 'genes.fasta.gz'
      --output_folder       Location for output
      --taxonomic_dmnd      Database used for taxonomic annotation (default: false)
      --ncbi_taxdump        Reference describing the NCBI Taxonomy
                            (default: ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
      --eggnog_dmnd         One of two databases used for functional annotation with eggNOG (default: false)
      --eggnog_db           One of two databases used for functional annotation with eggNOG (default: false)


    ####################################
    # Downloading Reference Databases: #
    ####################################

    --taxonomic_dmnd
        The DIAMOND database of reference protein sequences must be indexed using both
        (a) a set of sequences to search and (b) taxonomic annotations for each.
        Full instructions for creating this indexed database file can be found
        here: https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#makedb-options

        Example:
        diamond makedb \
            --in <proteins_fasta> \
            --db <output_dmnd> \
            --taxonmap prot.accession2taxid.FULL.gz \
            --taxonnodes nodes.dmp \
            --taxonnames names.dmp

    --eggnog_dmnd & --eggnog_db
        The eggNOG database for functional annotation can be most easily downloaded
        using the edicated utility provided along with the eggNOG-mapper utility.
        The only flag which needs to be set when running the download utility is the
        destination folder for the downloaded files.

        Example:
        download_eggnog_data.py --data_dir data/
    
    """.stripIndent()
}


// Show help message if the user specifies the --help flag at runtime
if (params.help || params.gene_fasta == false || params.output_folder == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Show help message if the user does not specify any annotations
if (params.taxonomic_dmnd == false && params.eggnog_dmnd == false && params.eggnog_db == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

workflow {

    main: 
    
    // Make sure we can find the input files
    if(file(params.gene_fasta).isEmpty()){
        log.info"""Cannot find input file ${params.gene_fasta}""".stripIndent()
        exit 0
    }

    // Run the annotation steps on the gene catalog
    Annotation_wf(
        file(params.gene_fasta)
    )



}