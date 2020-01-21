#!/usr/bin/env nextflow

// Script to download reference genomes from NCBI

// Using DSL-2
nextflow.preview.dsl=2

// Set default values for parameters
params.manifest = false
params.output = false
params.help = false

// Import a single function for concatenation which will be used multiple times
include concatenateFiles as concatGenomes_1 from "./modules/general"
include concatenateFiles as concatGenomes_2 from "./modules/general"
include concatenateFiles as concatCSV_1 from "./modules/general"
include concatenateFiles as concatCSV_2 from "./modules/general"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/geneshot/make_reference_genomes <ARGUMENTS>

    Input File:   
    
    The input --manifest file expected by this tool is in the format provided
    by the Genomes portal from NCBI. Specifically, the table is comma-delimited
    and contains a column for GenBank FTP folder as well as "#Organism Name"
    and Assembly.
    
    Required Arguments:
      --manifest            Manifest TSV from the NCBI Genomes portal
      --output              Folder to write output files

    Output Files:

    All output files will be written to the --output folder. This includes:

      - reference_genomes.fasta.gz: The complete set of references in FASTA format
      - reference_genomes.headers.csv.gz: Table linking contig headers with genome names

    Reference:

      - https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.manifest == false || params.output == false){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Make sure that --output ends with a trailing "/" character
if (!params.output.endsWith("/")){
    output_folder = params.output.concat("/")
} else {
    output_folder = params.output
}

// Download the genomes from the FTP server
process downloadGenome {
    container "quay.io/fhcrc-microbiome/wget@sha256:98b90e8bb8a171182505f1e255b0bd85cbbda68f08c08b4877c3fc48e63ac82f"
    label "io_limited"
    errorStrategy "retry"
    
    input:
    tuple val(org_name), val(assembly_name), val(ftp_prefix)

    output:
    file "${assembly_name.replaceAll(/"/, "")}.fasta.gz"
    file "${assembly_name.replaceAll(/"/, "")}.fasta.gz.csv.gz"

    """
# Break on any errors
set -e

# Download from the specified folder
ftp_prefix="${ftp_prefix.replaceAll(/"/, "")}"
ftp_suffix="\$(echo \$ftp_prefix | sed 's/.*\\///')"
genome_suffix=_genomic.fna.gz
url=\$ftp_prefix/\$ftp_suffix\$genome_suffix

# Set the name of the downloaded file
genome_fp="${assembly_name.replaceAll(/"/, "")}.fasta.gz"

echo Fetching \$url
wget \$url -O \$genome_fp || \
wget \$url -O \$genome_fp

# Make sure the resulting file is gzipped
gzip -t \$genome_fp

# Make sure that the file isn't empty
(( \$(gunzip -c \$genome_fp | wc -l) > 0 ))

# Strip everything after the first whitespace from the header
gunzip -c \$genome_fp | sed 's/ .*//' | gzip -c > TEMP && \
mv TEMP \$genome_fp

# Make a headerless CSV with the first column containing the
# headers from this FASTA, and the second containing a name
# for this organism which includes both the #Organism Name
# as well as the Assembly
genome_name=\$( echo ${org_name}_${assembly_name} | sed 's/[^A-Za-z0-9.]/_/g' )

gunzip -c \$genome_fp | grep '>' | tr -d '>' | while read header; do
    echo \$header,\$genome_name
done | gzip -c > \$genome_fp.csv.gz
    """

}

workflow {

    // Read four columns from the manifest CSV
    //  - #Organism Name
    //  - Assembly
    //  - GenBank FTP

    Channel.from(
        file(params.manifest)
    ).splitCsv(
        header: true
    ).map {
        r -> [r["#Organism Name"], r["Assembly"], r["GenBank FTP"]]
    }.set {
        genome_ch
    }

    // Download the genomes
    downloadGenome(
        genome_ch
    )

    // Concatenate the genomes in two steps
    concatGenomes_1(
        downloadGenome.out[0].toSortedList().flatten().collate(100)
    )
    concatGenomes_2(
        concatGenomes_1.out.collect()
    )
    // Concatenate the CSV header key in two steps
    concatCSV_1(
        downloadGenome.out[1].toSortedList().flatten().collate(100)
    )
    concatCSV_2(
        concatCSV_1.out.collect()
    )

    publish:

    concatGenomes_2.out to: "${output_folder}", mode: "copy"
    concatCSV_2.out to: "${output_folder}", mode: "copy"
}