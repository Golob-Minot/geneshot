// Processes to perform de novo assembly and annotate those assembled sequences

// Default parameters
params.output_prefix = "geneshot"

// Containers
container__assembler = "quay.io/biocontainers/megahit:1.2.9--h8b12597_0"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

include { diamondDB } from "./alignment" params(
    output_folder: params.output_folder
)
include {
    linclust as linclustRound1;
    linclust as linclustRound2;
    linclust as linclustRound3;
    linclust as linclustRound4;
    linclust as linclustRound5;
    diamondDedup as dedupRound1;
    diamondDedup as dedupRound2;
    diamondDedup as dedupRound3;
    diamondDedup as dedupRound4;
    diamondDedup as dedupRound5;
 } from "./mmseqs" params(
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

    // Calculate summary metrics for every assembled gene in each sample
    parseGeneAnnotations(
        prodigal.out[0]
    )

    // Combine genes by amino acid identity in five rounds
    // Each round will include both linclust- and DIAMOND-based deduplication
    // linclust provides fast symmetrical overlap search, while
    // DIAMOND performs slower, asymmetrical overlap search
    linclustRound1(
        prodigal.out[0].map{ r -> r[1]}.toSortedList().flatten().collate(4)
    )
    dedupRound1(
        linclustRound1.out
    )
    linclustRound2(
        dedupRound1.out.toSortedList().flatten().collate(4)
    )
    dedupRound2(
        linclustRound2.out
    )
    linclustRound3(
        dedupRound2.out.toSortedList().flatten().collate(4)
    )
    dedupRound3(
        linclustRound3.out
    )
    linclustRound4(
        dedupRound3.out.toSortedList().flatten().collate(4)
    )
    dedupRound4(
        linclustRound4.out
    )
    linclustRound5(
        dedupRound4.out.toSortedList()
    )
    dedupRound5(
        linclustRound5.out
    )

    // Make new, shorter names for each gene
    renameGenes(
        dedupRound5.out
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

    // Join the gene annotation tables with the gene assignments
    annotateAssemblies(
        parseGeneAnnotations.out.join(
            alignAlleles.out
        )
    )

    emit:
        gene_fasta = renameGenes.out
        allele_assembly_csv_list = annotateAssemblies.out.toSortedList()

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
cat OUTPUT/final.contigs.fa | sed 's/>/>${specimen}__GENE__/' | sed 's/ /__/g' | gzip -c > "${specimen}.contigs.fasta.gz"

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
        tuple val(specimen), file("${specimen}.faa.gz")
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
process parseGeneAnnotations {
    tag "Summarize every assembled gene"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'retry'
    
    input:
    tuple val(specimen), file(faa)
    
    output:
    tuple val(specimen), file("${specimen}.gene_annotations.csv.gz")

"""
#!/usr/bin/env python3
import gzip
import pandas as pd

# Function to parse a FASTA file with information encoded by Prodigal and Megahit in the header
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

    return pd.DataFrame(dat)

# Function to parse a single header
def parse_header(line):

    # The header line follows a format somewhat like this:
    # >Mock__15__GENE__k21_0__flag=1__multi=20.5919__len=48411_1 # 510 # 716 # -1 # ID=1_1;partial=00;start_type=ATG;rbs_motif=GGAG/GAGG;rbs_spacer=5-10bp;gc_cont=0.386
    #  ------------------------GENE_NAME------------------------   START STOP  STRAND
    #  SPECIMEN        CONTIG -----------NAME_DETAILS-----------                    -----------------------------------HEADER_DETAILS-----------------------------------

    # First let's just get the gene name
    gene_name, header = line[1:].rstrip("\\n").lstrip(">").split(" ", 1)

    # The "__GENE__" separator is used to preserve the specimen name within the gene name
    specimen, gene_remainder = gene_name.split("__GENE__", 1)

    # The contig name is the next field encoded in the gene name
    contig, gene_remainder = gene_remainder.split("__", 1)

    # The rest of the gene details are encoded in the gene_remainder (delimited by __)
    # as well as the final string in the header (delimited by ;)

    # We can parse the gene_remainder details with the __ delimiter and the =
    # We have to make sure to strip off the gene index number

    output_dat = dict([
        (field.split("=",1)[0], field.split("=",1)[1].split("_", 1)[0])
        for field in gene_remainder.split("__")
    ])

    # Now let's add in the details from the header
    start, stop, strand, header = header.strip(" ").strip("#").split(" # ")

    # Iterate over each of the elements provided
    for f in header.split(";"):
        if "=" in f:
            k, v = f.split("=", 1)
            output_dat[k] = v
    assert "gc_cont" in output_dat, (output_dat, line)

    # Add the other metrics to the output
    output_dat["gene_name"] = gene_name
    output_dat["start"] = int(start)
    output_dat["stop"] = int(stop)
    output_dat["strand"] = strand
    output_dat["specimen"] = specimen
    output_dat["contig"] = contig

    # Make some customizations to the data types
    del output_dat["ID"]
    for k, t in [
        ("strand", int), 
        ("len", int), 
        ("multi", float), 
        ("gc_cont", float)
    ]:
        assert k in output_dat, (output_dat, line)
        output_dat[k] = t(output_dat[k])

    return output_dat

# Parse and write out to a file
parse_prodigal_faa("${faa}").to_csv(
    "${specimen}.gene_annotations.csv.gz",
    index = None,
    compression = "gzip"
)


"""
}

// Combine the assembly table with the assignment of catalog genes
process annotateAssemblies {
    tag "Summarize every assembled gene"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'retry'

    publishDir "${params.output_folder}/assembly/${specimen}", mode: "copy"
    
    input:
    tuple val(specimen), file(assembly_csv), file(alignment_tsv)
    
    output:
    file "${specimen}.csv.gz"

"""
#!/usr/bin/env python3
import pandas as pd

print("Reading in assembly information for ${specimen}")
assembly_df = pd.read_csv("${assembly_csv}")
print("Read in %d assembled genes" % assembly_df.shape[0])

# Read in the table linking each allele to the gene catalog
gene_alignments = pd.read_csv(
    "${alignment_tsv}",
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

# Add the gene name to the gene annotation table
assembly_df["catalog_gene"] = assembly_df["gene_name"].apply(
    gene_alignments.groupby(
        "allele"
    ).head(
        1
    ).set_index(
        "allele"
    )[
        "gene"
    ].get
)
print("%d / %d genes have an annotation in the gene catalog" % (assembly_df["catalog_gene"].dropna().shape[0], assembly_df.shape[0]))

# Save and exit
assembly_df.to_csv(
    "${specimen}.csv.gz",
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
    --temp_dir TEMP/

gzip genes.emapper.annotations
    
    """

}


// Assign a new, shorter name to a set of genes
process renameGenes {
    tag "Make concise unique gene names"
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'io_limited'
    errorStrategy 'retry'
    publishDir "${params.output_folder}/ref/", mode: "copy"

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

used_strings = set([])

with gzip.open("genes.fasta.gz", "wt") as fo:
    with gzip.open("input.genes.fasta.gz", "rt") as fi:
        for header, seq in SimpleFastaParser(fi):
            new_string = random_string()
            while new_string in used_strings:
                new_string = random_string()
            used_strings.add(new_string)
            fo.write(">gene_%s_%daa\\n%s\\n" % (new_string, len(seq), seq))

"""
}

// Use alignment to figure out which assembled allele was grouped into which gene
process alignAlleles {
    tag "Match alleles to gene centroids"
    container "quay.io/fhcrc-microbiome/famli:v1.5"
    label 'mem_medium'
    errorStrategy 'retry'
    
    input:
    tuple val(specimen), file(alleles_fasta)
    file refdb
    
    output:
    tuple val(specimen), file("${specimen}.gene_alignments.tsv.gz") optional true

    """
    set -e

    echo "Aligning ${alleles_fasta}"

    # Check to see if there are any reads
    if (( \$( gunzip -c ${alleles_fasta} | wc -l ) <= 1 )); then
    
        echo "No alleles found in ${alleles_fasta}, skipping"
    
    else

        # Make the output filepath
        fo="${specimen}.gene_alignments.tsv.gz"
        echo "Writing out to \$fo"

        diamond \
        blastp \
        --query ${alleles_fasta} \
        --out \$fo \
        --threads ${task.cpus} \
        --db ${refdb} \
        --outfmt 6 qseqid sseqid pident length qlen slen \
        --query-cover ${params.min_coverage} \
        --id ${params.min_identity} \
        --top 0 \
        --block-size ${task.memory.toMega() / (1024 * 6)} \
        --compress 1 \
        --unal 0

    fi
    """

}
