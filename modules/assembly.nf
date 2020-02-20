// Processes to perform de novo assembly and annotate those assembled sequences

// Containers
container__assembler = "quay.io/biocontainers/megahit:1.2.9--h8b12597_0"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas@sha256:b57953e513f1f797522f88fa6afca187cdd190ca90181fa91846caa66bdeb5ed"

include diamondDB from "./alignment" params(
    output_folder: params.output_folder
)
include linclust as linclustRound1 from "./mmseqs" params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include linclust as linclustRound2 from "./mmseqs" params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include linclust as linclustRound3 from "./mmseqs" params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include linclust as linclustRound4 from "./mmseqs" params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)
include linclust as linclustRound5 from "./mmseqs" params(
    min_identity: params.min_identity,
    min_coverage: params.min_coverage
)

workflow assembly_wf {
    take:
        combined_reads_ch

    main:

    // Perform de novo assembly
    assembly(
        combined_reads_ch
    )

    // Annotate those contigs with Prokka
    prodigal(
        assembly.out
    )

    // Calculate summary metrics for every assembled gene in batches of 10 samples
    geneAssemblyMetrics(
        prodigal.out[0].toSortedList().flatten().collate(10)
    )

    // Join together each batch of assembly summaries
    joinAssemblyMetrics(
        geneAssemblyMetrics.out.toSortedList()
    )

    // Combine genes by amino acid identity in five rounds
    linclustRound1(
        prodigal.out[0].toSortedList().flatten().collate(4)
    )
    linclustRound2(
        linclustRound1.out.toSortedList().flatten().collate(4)
    )
    linclustRound3(
        linclustRound2.out.toSortedList().flatten().collate(4)
    )
    linclustRound4(
        linclustRound3.out.toSortedList().flatten().collate(4)
    )
    linclustRound5(
        linclustRound4.out.toSortedList()
    )

    // Make new, shorter names for each gene
    renameGenes(
        linclustRound5.out
    )

    // Index the assembled alleles for alignment
    diamondDB(
        renameGenes.out
    )

    // Align the assembled alleles against the gene centroids
    alignAlleles(
        prodigal.out[0],
        diamondDB.out
    )

    // Make a single table with all of the allele - gene pairings
    makeAlleleTable(
        alignAlleles.out.toSortedList()
    )

    emit:
        gene_fasta = renameGenes.out
        allele_gene_tsv = makeAlleleTable.out
        allele_assembly_csv = joinAssemblyMetrics.out

}

// Break out the annotation into a separate workflow
workflow annotation_wf {
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
        join_eggnog(
            eggnog_tsv
        )
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


// De novo assembly
process assembly {
    tag "De novo metagenomic assembly"
    container "${container__assembler}"
    label 'mem_veryhigh'
    errorStrategy "retry"

    publishDir "${params.output_folder}/assembly/${specimen}", mode: "copy"

    input:
        tuple specimen, file(R1), file(R2)
    
    output:
        tuple specimen, file("${specimen}.contigs.fasta.gz"), file("${specimen}.megahit.log")
    
"""
set -e 

date
echo -e "Running Megahit\\n"

megahit \
    -1 ${R1} -2 ${R2} \
    -o OUTPUT \
    -t ${task.cpus}

date
echo -e "\\nMaking sure output files are not empty\\n"
[[ \$(cat OUTPUT/final.contigs.fa | wc -l) > 0 ]]

date
echo -e "\\nRenaming output files\\n"

mv OUTPUT/log "${specimen}.megahit.log"

# Add the specimen name to the contig name
cat OUTPUT/final.contigs.fa | sed 's/>/>${specimen}__/' | sed 's/ /__/g' | gzip -c > "${specimen}.contigs.fasta.gz"

date
echo -e "\\nDone\\n"
"""
}

// Annotation of coding sequences with prodigal

process prodigal {
    tag "Identify protein-coding genes"
    container 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
    label 'io_limited'
    errorStrategy "retry"
    publishDir "${params.output_folder}/assembly/${specimen}/", mode: "copy"

    input:
        tuple val(specimen), file(contigs), file(spades_log)
    
    output:
        file "${specimen}.faa.gz"
        file "${specimen}.gff.gz"
    
"""
set -e 

gunzip -c ${contigs} > ${specimen}.contigs.fasta

prodigal \
    -a ${specimen}.faa \
    -i  ${specimen}.contigs.fasta \
    -f gff \
    -o ${specimen}.gff \
    -p meta

gzip ${specimen}.gff
gzip ${specimen}.faa

"""
}


// Summarize the depth of sequencing and GC content for every assembled gene
process geneAssemblyMetrics {
    tag "Summarize every assembled gene"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'retry'
    
    input:
    file faa_list
    
    output:
    file "allele.assembly.metrics.csv.gz"

"""
#!/usr/bin/env python3
import gzip
import pandas as pd

# Parse the list of FASTA files to read in
faa_list = "${faa_list}".split(" ")

# Function to parse a FASTA file with information encoded by Prodigal and metaSPAdes in the header
def parse_prodigal_faa(fp):

    # Keep track of all of the information for each assembled gene
    dat = []

    # Open a connection to the file
    with gzip.open(fp, "rt") as f:

        # Iterate over every line
        for line in f:

            # Skip the non-header lines
            if line.startswith(">") is False:
                continue

            # Add the information for this header
            dat.append(parse_header(line))

    return dat

# Function to parse a single header
def parse_header(line):

    # Keep track of the gene name
    # The next two fields after the gene name are the start and the stop positions
    gene_name, start, stop, strand, details = line.rstrip("\\n").split(" # ")

    # The first field in the gene name is the specimen
    specimen, gene_name = gene_name[1:].split("__", 1)

    # There is more good information to get from the gene name
    output_dat = dict([
        (field.split("=",1)[0], field.split("=",1)[1])
        for field in gene_name.split("__")
        if "=" in field
    ])

    # The final piece we want to get is the GC content, from the `details`
    gc = None
    for f in details.split(";"):
        if f.startswith("gc_cont="):
            gc = f[len("gc_cont="):]
    assert gc is not None

    # Add the other metrics to the output
    output_dat["gene_name"] = gene_name
    output_dat["start"] = int(start)
    output_dat["stop"] = int(stop)
    output_dat["strand"] = strand
    output_dat["details"] = details
    output_dat["specimen"] = specimen
    output_dat["gc"] = float(gc)

    return output_dat

# Read in all of the files
df = pd.DataFrame([
    gene_details
    for fp in faa_list
    for gene_details in parse_prodigal_faa(fp)
])

# Write out to a file
df.to_csv(
    "allele.assembly.metrics.csv.gz",
    index = None,
    compression = "gzip"
)


"""
}

// Join together the assembly summaries from each batch of 10 samples processed earlier
process joinAssemblyMetrics {
    tag "Summarize every assembled gene"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'retry'
    publishDir "${params.output_folder}/assembly/", mode: "copy"
    
    input:
    file "allele.assembly.metrics.*.csv.gz"
    
    output:
    file "allele.assembly.metrics.csv.gz"

"""
#!/usr/bin/env python3
import gzip
import os
import pandas as pd

# Parse the list of CSV files to read in
csv_list = [fp for fp in os.listdir(".") if fp.startswith("allele.assembly.metrics.") and fp.endswith(".csv.gz")]

print("Preparing to read in %d CSV files" % len(csv_list))

# Read in all of the CSV files
df = pd.concat([
    pd.read_csv(fp, sep=",", compression="gzip")
    for fp in csv_list
])

# Write out to a file
df.to_csv(
    "allele.assembly.metrics.csv.gz",
    index = None,
    compression = "gzip"
)


"""
}


process shard_genes {
    tag "Split the gene catalog into smaller shards"
    container "ubuntu:18.04"
    label 'mem_medium'
    errorStrategy 'retry'
    
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


process diamond_tax {
    tag "Annotate genes by taxonomy"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
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
    container "quay.io/biocontainers/eggnog-mapper:2.0.1--py_1"
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
    --output genes \
    -m "diamond" \
    --cpu ${task.cpus} \
    --data_dir data/ \
    --scratch_dir SCRATCH/ \
    --temp_dir TEMP/ \

gzip genes.emapper.annotations
    
    """

}

process join_eggnog {
    tag "Concatenate eggNOG annotation files"
    container "ubuntu:18.04"
    label 'mem_medium'
    publishDir "${params.output_folder}/annot/", mode: "copy"

    input:
    file "genes.emapper.annotations.*.gz"
    
    output:
    file "genes.emapper.annotations.gz"

    
"""
set -e

for fp in genes.emapper.annotations.*.gz; do

    cat \$fp
    rm \$fp

done > genes.emapper.annotations.gz

"""

}

// Assign a new, shorter name to a set of genes
process renameGenes {
    tag "Make concise unique gene names"
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'io_limited'
    errorStrategy 'retry'

    input:
    file "input.genes.fasta.gz"

    output:
    file "genes.fasta.gz"

"""
#!/usr/bin/env python3

from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import uuid

def random_string(n=8):
    return str(uuid.uuid4())[:n]

with gzip.open("genes.fasta.gz", "wt") as fo:
    with gzip.open("input.genes.fasta.gz", "rt") as fi:
        for header, seq in SimpleFastaParser(fi):
            fo.write(">gene_%s_%daa\\n%s\\n" % (random_string(), len(seq), seq))

"""
}

// Use alignment to figure out which assembled allele was grouped into which gene
process alignAlleles {
    tag "Match alleles to gene centroids"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_medium'
    errorStrategy 'retry'
    
    input:
    file alleles_fasta
    file refdb
    
    output:
    file "${alleles_fasta}.tsv.gz"

    """
    set -e

    diamond \
      blastp \
      --query ${alleles_fasta} \
      --out ${alleles_fasta}.tsv.gz \
      --threads ${task.cpus} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length qlen slen \
      --query-cover ${params.min_coverage} \
      --id ${params.min_identity} \
      --top 0 \
      --block-size ${task.memory.toMega() / (1024 * 6)} \
      --compress 1 \
      --unal 0
    """

}

// Make a single table linking every assembledallele to every gene centroid
process makeAlleleTable {
    tag "Combine allele tables across samples"
    container "${container__pandas}"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    file tsv_list
    
    output:
    file "alleles.genes.tsv.gz"

"""
#!/usr/bin/env python3

import pandas as pd
import os

# Function to read in the allele tables
def read_allele_table(fp):

    assert os.path.exists(fp)

    print("Reading in %s" % fp)

    return pd.read_csv(
        fp,
        sep = "\\t",
        header = None,
        names = [
            "allele",
            "gene",
            "pct_iden",
            "alignment_len",
            "allele_len",
            "gene_len"
        ],
        compression = "gzip"
    )

# Make a single table and write it out
pd.concat([
    read_allele_table(fp)
    for fp in "${tsv_list}".split(" ")
]).to_csv(
    "alleles.genes.tsv.gz",
    sep = "\\t",
    compression = "gzip",
    index = None
)
"""

}
