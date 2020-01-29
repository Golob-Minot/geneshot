// Processes to perform de novo assembly and annotate those assembled sequences

// Containers
container__spades = "quay.io/biocontainers/spades:3.14.0--h2d02072_0"

include makeDiamondDB from "./alignment"

workflow assembly_wf {
    get:
        combined_reads_ch

    main:

    // Assemble samples with metaSPAdes
    metaspadesAssembly(
        combined_reads_ch
    )

    // Annotate those contigs with Prokka
    prodigalAnnotate(
        metaspadesAssembly.out
    )

    // Calculate summary metrics for every assembled gene
    geneAssemblyMetrics(
        prodigalAnnotate.out[0].collect()
    )

    // Combine the gene sequences across all samples
    combineCDS(
        prodigalAnnotate.out[0].collect()
    )

    // Combine genes by amino acid identity
    clusterCDS(
        combineCDS.out
    )

    emit:
        gene_fasta = clusterCDS.out[0]
        allele_gene_tsv = clusterCDS.out[1]
        allele_assembly_csv = geneAssemblyMetrics.out

}

// Break out the annotation into a separate workflow
workflow annotation_wf {
    get:
    gene_fasta

    main:

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
        eggnog_annotation(
            gene_fasta,
            file(params.eggnog_db),
            file(params.eggnog_dmnd)
        )
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
        taxonomic_annotation(
            gene_fasta,
            file(params.taxonomic_dmnd)
        )
    }

    // Determine whether or not to run the genome alignment based
    // on --noannot, --ref_genome_fasta, and --ref_genome_csv
    run_genome_alignment = false
    if ( params.noannot == false ) {
        if ( params.ref_genome_fasta && params.ref_genome_csv ) {
            if ( !file(params.ref_genome_fasta).isEmpty() && !file(params.ref_genome_csv).isEmpty() ){
                run_genome_alignment = true
            }
        }
    }

    // Annotate the clustered genes by alignment against whole reference genomes
    if ( run_genome_alignment ) {

        // Make a DIAMOND database for the provided genes
        makeDiamondDB(
            gene_fasta
        )

        alignGenomes(
            makeDiamondDB.out,
            file(params.ref_genome_fasta)
        )
    }

}


// Assembly with metaspades
process metaspadesAssembly {
    tag "De novo assembly with metaSPAdes"
    container "${container__spades}"
    label 'mem_veryhigh'
    errorStrategy "retry"

    publishDir "${params.output_folder}/assembly/${specimen}", mode: "copy"

    input:
        tuple specimen, file(R1), file(R2)
    
    output:
        tuple specimen, file("${specimen}.contigs.fasta.gz"), file("${specimen}.scaffolds.fasta.gz"), file("${specimen}.metaspades.log")
    
"""
set -e 

metaspades.py \
    --meta \
    --phred-offset ${params.phred_offset} \
    -1 ${R1} -2 ${R2} \
    -o . \
    -t ${task.cpus} -m ${task.memory.toMega() / 1024} | 
    tee -a ${specimen}.metaspades.log

echo "\\n\\nDone\\n\\n"

# Print the spades log to help with troubleshooting in the future
echo "Log file:\\n\\n"
cat spades.log

# Make sure there were no errors in the log
(( \$(cat ${specimen}.metaspades.log | grep -c "== Error ==" ) == 0 ))

# Make sure that the contig files aren't empty
[[ -s contigs.fasta ]] && (( \$(cat contigs.fasta | wc -l) > 1))
[[ -s scaffolds.fasta ]] && (( \$(cat scaffolds.fasta | wc -l) > 1))

# Add the specimen name to the contig name
cat contigs.fasta | sed 's/>/>${specimen}_/' | gzip -c > ${specimen}.contigs.fasta.gz
cat scaffolds.fasta | sed 's/>/>${specimen}_/' | gzip -c > ${specimen}.scaffolds.fasta.gz
"""
}

// Annotation of coding sequences with prodigal

process prodigalAnnotate {
    tag "Identify protein-coding genes"
    container 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
    label 'io_limited'
    errorStrategy "retry"
    publishDir "${params.output_folder}/assembly/${specimen}/", mode: "copy"

    input:
        tuple val(specimen), file(contigs), file(scaffolds), file(spades_log)
    
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
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    label 'io_limited'
    errorStrategy 'retry'
    publishDir "${params.output_folder}/assembly/", mode: "copy"
    
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
    gene_name, start, stop, strand, details = line.split(" # ")

    # The first field in the gene name is the specimen
    specimen, gene_name = gene_name[1:].split("_NODE_", 1)

    # There is more good information to get from the gene name
    contig_num, _, length, _, depth, gene_num = gene_name.split("_")

    # The final piece we want to get is the GC content, from the `details`
    gc = None
    for f in details.split(";"):
        if f.startswith("gc_cont="):
            gc = f[len("gc_cont="):]
    assert gc is not None

    return dict([
        ("gene_name", gene_name),
        ("start", int(start)),
        ("stop", int(stop)),
        ("strand", strand),
        ("details", details),
        ("specimen", specimen),
        ("contig_num", contig_num),
        ("contig_length", float(length)),
        ("contig_depth", float(depth)),
        ("gene_num", gene_num),
        ("gc", float(gc))
    ])

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


process combineCDS {
    tag "Combine gene sequences"
    container "ubuntu:16.04"
    label 'io_limited'
    errorStrategy 'retry'
    maxRetries 10
    
    input:
    file "assembly.*.fasta.gz"
    
    output:
    file "all_CDS.fasta.gz"

    """
#!/bin/bash
set -e
cat *fasta.gz > all_CDS.fasta.gz
gzip -t all_CDS.fasta.gz
    """
}


process clusterCDS {
    tag "Cluster genes with similar sequences"
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    publishDir "${params.output_folder}/ref/", mode: "copy"
    
    input:
    file all_cds
    
    output:
    file "genes.fasta.gz"
    file "genes.alleles.tsv.gz"
    
    """
#!/bin/bash
set -e
# Make the MMSeqs2 database
mmseqs createdb ${all_cds} db
# Cluster the protein sequences
mmseqs linclust db cluster_db ./ \
    --min-seq-id ${params.min_identity / 100} \
    --max-seqs 100000 \
    -c ${params.min_coverage / 100}
# Make TSV output for clustering
mmseqs createtsv db db cluster_db genes.alleles.tsv
# Get the representative sequences
mmseqs result2repseq db cluster_db genes
mmseqs result2flat db db genes genes.fasta --use-fasta-header
gzip genes.alleles.tsv
gzip genes.fasta
    """
}


process taxonomic_annotation {
    tag "Annotate genes by taxonomy"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'
    publishDir "${params.output_folder}/annot/", mode: "copy"

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


process eggnog_annotation {
    tag "Annotate genes by predicted function"
    container "quay.io/biocontainers/eggnog-mapper:2.0.1--py_1"
    label 'mem_veryhigh'
    publishDir "${params.output_folder}/annot/", mode: "copy"
    
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

process alignGenomes {
    tag "Align genes against reference genomes"
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label "mem_veryhigh"
    // errorStrategy "retry"
    
    input:
    file genes_dmnd
    file fasta_gz

    output:
    file "reference_genomes.aln.gz"

    """
set -e

# Make sure all files are present
echo "Making sure ${genes_dmnd} is present"
[[ -s ${genes_dmnd} ]]
echo "Making sure ${fasta_gz} is present"
[[ -s ${fasta_gz} ]]

echo "Processing ${fasta_gz} and ${genes_dmnd}"

diamond \
    blastx \
    --query ${fasta_gz} \
    --out reference_genomes.aln.gz \
    --threads ${task.cpus} \
    --db ${genes_dmnd} \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
    --range-cover ${params.min_coverage} \
    --subject-cover ${params.min_coverage} \
    --range-culling \
    --id ${params.min_identity} \
    --top 1 \
    --block-size ${task.memory.toMega() / (1024 * 6)} \
    --query-gencode ${params.gencode} \
    -F 15 \
    --unal 0 \
    --compress 1

echo "Done"
    """

}
