// Processes to perform de novo assembly and annotate those assembled sequences
nextflow.enable.dsl=2
// Default parameters
params.output_prefix = "geneshot"
// Assembly options
params.gene_fasta = false
params.phred_offset = 33 // spades
params.min_identity = 50 // linclust and reference genome alignment
params.min_coverage = 80 // linclust and reference genome alignment
// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.help = false
params.output = './results/'
params.manifest = null

// Make sure that --output ends with trailing "/" characters
if (!params.output.endsWith("/")){
    params.output_folder = params.output.concat("/")
} else {
    params.output_folder = params.output
}

// Containers
container__assembler = "quay.io/biocontainers/megahit:1.2.9--h8b12597_0"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"
container__prodigal = 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
container__fastatools = "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.2"

include { DiamondDB } from "./alignment" params(
    output_folder: params.output_folder
)
include { MMSeqs2_Cluster as MMSeqs2_Cluster_100 } from "./mmseqs" params(
    min_identity: 100,
    min_coverage: 80
)
include { MMSeqs2_Cluster as MMSeqs2_Cluster_90 } from "./mmseqs" params(
    min_identity: 90,
    min_coverage: 80
)
include { MMSeqs2_Cluster as MMSeqs2_Cluster_50 } from "./mmseqs" params(
    min_identity: 50,
    min_coverage: 80
)


workflow Genecatalog_wf {
    take:
        combined_reads_ch

    main:

    // Perform de novo assembly
    Assembly(
        combined_reads_ch
    )

    // Annotate those contigs with Prokka
    Prodigal(
        Assembly.out
    )

    // Calculate summary metrics for every assembled gene in each sample
    ParseGeneAnnotations(
        Prodigal.out[0].map{ r -> [r[0], r[1]]}
    )

    // Make a dereplicated allele catalog 
    DereplicateAlleles(
        ParseGeneAnnotations.out.collect{ it[0] }, // specimens
        ParseGeneAnnotations.out.collect{ it[1] }, // specimen allele FAA
        ParseGeneAnnotations.out.collect{ it[2] }, // specimen allele info CSV
    )

    // Cluster at 80 / 100 c/i
    MMSeqs2_Cluster_100(
        DereplicateAlleles.out.alleles
    )

    MMSeqs2_Cluster_90(
        MMSeqs2_Cluster_100.out.seqs
    )

    MMSeqs2_Cluster_50(
        MMSeqs2_Cluster_90.out.seqs
    )
    // Summarize!

    SummarizeAllelesAndClusters(
        DereplicateAlleles.out.allele_info,
        MMSeqs2_Cluster_100.out.clusters,
        MMSeqs2_Cluster_90.out.clusters,
        MMSeqs2_Cluster_50.out.clusters,
        MMSeqs2_Cluster_100.out.seqs,
        MMSeqs2_Cluster_90.out.seqs,
        MMSeqs2_Cluster_50.out.seqs,        
    )

    //emit:
    //    gene_fasta = renameGenes.out
    //    allele_assembly_csv_list = AnnotateAssemblies.out.map{ it[1] }.toSortedList()

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
process Assembly {
    tag "De novo metagenomic assembly"
    container "${container__assembler}"
    label 'mem_veryhigh'
    errorStrategy "ignore"

    publishDir "${params.output_folder}/genecatalog/${specimen}", mode: "copy"

    input:
        tuple val(specimen), file(R1), file(R2)
    
    output:
        tuple val(specimen), file("${specimen}.contigs.fasta.gz"), file("${specimen}.megahit.log")
    
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

process Prodigal {
    tag "Identify protein-coding genes"
    container "${container__prodigal}"
    label 'mem_medium'
    errorStrategy "finish"
    publishDir "${params.output_folder}/genecatalog/${specimen}/", mode: "copy"

    input:
        tuple val(specimen), file(contigs), file(spades_log)
    
    output:
        tuple val(specimen), file("${specimen}.faa.gz"), file("${specimen}.fna.gz")
        file "${specimen}.gff.gz"
    
"""
set -e 

gunzip -c ${contigs} > ${specimen}.contigs.fasta

prodigal \
    -a ${specimen}.faa \
    -d ${specimen}.fna \
    -i  ${specimen}.contigs.fasta \
    -f gff \
    -o ${specimen}.gff \
    -p meta

gzip ${specimen}.gff
gzip ${specimen}.faa
gzip ${specimen}.fna

"""
}


// Dereplicate and output alleles
process DereplicateAlleles {
    tag "Dereplicate and output *alleles*"
    container "${container__fastatools}"
    label 'mem_medium'
    errorStrategy 'finish'
    
    input:
    val(specimens)
    file(faas)
    file(specimen_allele_csvs)
    
    output:
    path "allele_info.csv.gz", emit: allele_info
    path "alleles.faa.gz", emit: alleles

    publishDir "${params.output_folder}/genecatalog/", mode: "copy"

"""
#!/usr/bin/env python3
import gzip
from collections import defaultdict
import fastalite
import csv

# First load in the per_specimen_allele -> contig
sp_allele_csvs = "${specimen_allele_csvs.join(';;')}".split(';;')
specimen_allele_info = {}
for sac in sp_allele_csvs:
    sac_r = csv.DictReader(gzip.open(sac, 'rt'))
    specimen_allele_info.update({
        r['gene_name']: {
            'specimen_allele': r['gene_name'],
            'specimen': r['specimen'],
            'contig': r['contig']
        } 
        for r in sac_r
    })
specimens = "${specimens.join(';;')}".split(';;')
allele_faas = "${faas.join(';;')}".split(';;')

centroid_allele = defaultdict(set)

for faa in allele_faas:
    for sr in fastalite.fastalite(gzip.open(faa, 'rt')):
        centroid_allele[sr.seq].add(sr.id)

# Great, now rename and output alleles into fasta format
with gzip.open("alleles.faa.gz", 'wt') as allele_h:
    for centroid_n, (centroid_seq, seq_ids) in enumerate(centroid_allele.items()):
        centroid_id = "allele___{:010X}".format(centroid_n+1)
        for sID in seq_ids:
            specimen_allele_info[sID]['allele'] = centroid_id
        allele_h.write(">{}\\n{}\\n".format(
            centroid_id,
            centroid_seq
        ))

with gzip.open("allele_info.csv.gz", 'wt') as allele_info_h:
    allele_w = csv.DictWriter(allele_info_h, fieldnames=['allele', 'specimen_allele', 'specimen', 'contig'])
    allele_w.writeheader()
    allele_w.writerows(specimen_allele_info.values())
"""
}

// Combine the assembly table with the assignment of catalog genes
process AnnotateAssemblies {
    tag "Summarize every assembled gene"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'finish'

    publishDir "${params.output_folder}/genecatalog/${specimen}", mode: "copy"
    
    input:
    tuple val(specimen), file(faa), file(assembly_csv), file(alignment_tsv)
    
    output:
    tuple val(specimen), path("${specimen}.csv.gz"), path(faa)

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

// Summarize the depth of sequencing and GC content for every assembled gene
process ParseGeneAnnotations {
    tag "Summarize every assembled gene"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'finish'
    
    input:
    tuple val(specimen), file(faa)
    
    output:
    tuple val(specimen), file(faa), file("${specimen}.gene_annotations.csv.gz")

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

process SummarizeAllelesAndClusters {
    tag "Summarize all of the avaliable alleles and output"
    container "${container__pandas}"
    label 'mem_medium'
    errorStrategy 'finish'
    publishDir "${params.output_folder}/genecatalog/", mode: "copy"
    
    input:
        path Alleles_csv
        path C100_tsv
        path C90_tsv
        path C50_tsv
        path Centroids_C100
        path Centroids_C90
        path Centroids_C50
    
    output:
        path('Allele_Cluster_info.csv.gz')
        path Centroids_C100
        path Centroids_C90
        path Centroids_C50        

"""
#!/usr/bin/env python3
import pandas as pd
a_i = pd.read_csv('${Alleles_csv}')
A_c100 = {
    r.allele: r.C100
    for i, r in 
    pd.read_csv('${C100_tsv}', delimiter='\\t', names=['C100', 'allele']).iterrows()
}
a_i['C100'] = a_i.allele.apply(A_c100.get)
c100_c90 = {
    r.C100: r.C90
    for i, r in 
    pd.read_csv('${C90_tsv}', delimiter='\\t', names=['C90', 'C100']).iterrows()
}
a_i['C90'] = a_i.C100.apply(c100_c90.get)
c90_c50 = {
    r.C90: r.C50
    for i, r in 
    pd.read_csv('${C50_tsv}', delimiter='\\t', names=['C50', 'C90']).iterrows()
}
a_i['C50'] = a_i.C90.apply(c90_c50.get)
a_i.to_csv(
    "Allele_Cluster_info.csv.gz",
    index=None
)


"""
}

// ---



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
    errorStrategy 'finish'
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
process AlignAlleles {
    tag "Match alleles to gene centroids"
    container "quay.io/fhcrc-microbiome/famli:v1.5"
    label 'mem_medium'
    errorStrategy 'finish'
    
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
        blastx \
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

//
// Steps to run gene-catalog independently.
//


// imports
include { Read_manifest } from './general'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run Golob-Minot/geneshot/modules/genecatalog.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    Options:
      --output              Folder to place analysis outputs (default ./results)
      --output_prefix       Text used as a prefix for summary HDF5 output files (default: geneshot)
      -w                    Working directory. Defaults to `./work`

    For Assembly:
      --phred_offset        for spades. Default 33.
      --min_identity        Amino acid identity cutoff used to combine similar genes (default: 50)
      --min_coverage        Length cutoff used to combine similar genes (default: 80) (linclust)

    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This CANNOT be repeated. One row per specimen
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      These MUST be preprocessed for this workflow.
    """.stripIndent()
}

workflow {
    main:

 
    // Show help message if the user specifies the --help flag at runtime
    if (params.help || params.manifest == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }
    // Read and validate manifest
    manifest_file = Channel.from(file(params.manifest))
    manifest = Read_manifest(manifest_file).valid_paired


    Genecatalog_wf(
        manifest.map{ [it.specimen, file(it.R1), file(it.R2)] }
  )

}
