#!/usr/bin/env nextflow

// Processes to perform de novo assembly and annotate those assembled sequences
nextflow.enable.dsl=2
// Default parameters
// Assembly options
params.gene_fasta = false
params.phred_offset = 33 // spades

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
container__anndata = "golob/python-anndata:0.9.2"
container__prodigal = 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
container__fastatools = "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.2"

include { DiamondDB as DiamondIndex_Alleles } from "./quantify" addParams(
    output_folder: params.output_folder
)

workflow Allele_catalog {
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

    DiamondIndex_Alleles(
        DereplicateAlleles.out.alleles
    )

    emit:
        alleles = DereplicateAlleles.out.alleles
        allele_info = DereplicateAlleles.out.allele_info
        alleles_dmdb = DiamondIndex_Alleles.out
    // */
}



// De novo assembly
process Assembly {
    tag "De novo metagenomic assembly"
    container "${container__assembler}"
    label 'mem_veryhigh'
    errorStrategy "ignore"

    publishDir "${params.output_folder}/contigs/${specimen}", mode: "copy"

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
    publishDir "${params.output_folder}/contigs/${specimen}/", mode: "copy"

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

    publishDir "${params.output_folder}/alleles/", mode: "copy"

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
            'contig': r['contig'],
            'strand': r['strand'],
            'start': r['start'],
            'stop': r['stop'],
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
    allele_w = csv.DictWriter(allele_info_h, fieldnames=['allele', 'specimen_allele', 'specimen', 'contig', 'strand', 'start', 'stop'])
    allele_w.writeheader()
    allele_w.writerows(specimen_allele_info.values())
"""
}

// Combine the assembly table with the assignment of catalog genes
process AnnotateAssemblies {
    tag "Summarize every assembled gene"
    container "${container__anndata}"
    label 'mem_medium'
    errorStrategy 'finish'

    publishDir "${params.output_folder}/alleles/${specimen}", mode: "copy"
    
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
    container "${container__anndata}"
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


//
// Steps to run allele-catalog independently.
//


// imports
include { Read_manifest } from './general'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run modules/allele_catalog.nf <ARGUMENTS>
    
    From a set of raw fastq reads, assemble contigs, extract protein coding sequences
    and dereplicate into an allele_catalog.
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    Options:
      --output              Folder to place analysis outputs (default ./results)
      -w                    Working directory. Defaults to `./work`

    For Assembly:
      --phred_offset        for spades. Default 33.

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


    Allele_catalog(
        manifest.map{ [it.specimen, file(it.R1), file(it.R2)] }
  )

}
