#!/usr/bin/env nextflow

// Script to download reference genomes from NCBI

// Using DSL-2
nextflow.preview.dsl=2

// Set default values for parameters
params.manifest = false
params.dbname = "reference_genomes"
params.output = false
params.help = false
params.batchsize = 100 // Number of genomes to download by a single worker

// Import a single function for concatenation which will be used multiple times
include concatenateFiles as concatGenomes_1 from "./modules/general"
include concatenateFiles as concatGenomes_2 from "./modules/general"
include concatenateFiles as concatCSV_1 from "./modules/general"
include concatenateFiles as concatCSV_2 from "./modules/general"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/geneshot/make_reference_genomes.nf <ARGUMENTS>

    Input File:   
    
    The input --manifest file expected by this tool is in the format provided
    by the Genomes portal from NCBI. Specifically, the table is comma-delimited
    and contains a column for GenBank FTP folder as well as "#Organism Name"
    and Assembly.
    
    Required Arguments:
      --manifest            Manifest TSV from the NCBI Genomes portal
      --output              Folder to write output files
      --dbname              Prefix for the output files (default: reference_genomes)

    Optional Arguments:
      --batchsize           Number of genomes to download by a single worker (default: 100)

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
// This process will download multiple genomes
// The set of genomes to download is contained in ${input_string_list}
// which contains a single string for each genome, formatted
// with tabs separating each field. Each string is separated by a newline.
process downloadGenome {
    container "quay.io/fhcrc-microbiome/wget@sha256:98b90e8bb8a171182505f1e255b0bd85cbbda68f08c08b4877c3fc48e63ac82f"
    label "io_limited"
    errorStrategy "retry"
    
    input:
    val genome_input_list

    output:
    file "ALL_GENOMES.fasta.gz"
    file "ALL_GENOMES.csv.gz"

    """
# Break on any errors
set -e

# Iterate over each genome
echo \"\"\"${genome_input_list}\"\"\" | while read org_name assembly_name ftp_prefix; do

    echo "Organism: \$org_name"
    echo "Assembly: \$assembly_name"
    echo "FTP Prefix: \$ftp_prefix"

    # Download from the specified folder
    ftp_suffix="\$(echo \$ftp_prefix | sed 's/.*\\///')"
    genome_suffix=_genomic.fna.gz
    url=\$ftp_prefix/\$ftp_suffix\$genome_suffix

    # Get the name of the genome
    genome_name=\$( echo "\${org_name}_\${assembly_name}" | sed 's/[^A-Za-z0-9.]/_/g' )

    # Set the name of the downloaded file
    genome_fp="\$genome_name.fasta.gz"

    echo Fetching \$url
    wget \$url -O \$genome_fp || \
    wget \$url -O \$genome_fp

    # Make sure the resulting file is gzipped
    gzip -t \$genome_fp

    # Make sure that the file isn't empty
    (( \$(gunzip -c \$genome_fp | wc -l) > 0 ))

    # Strip everything after the first whitespace from the header
    # Concatenate this genome to the set of ALL_GENOMES.fasta.gz
    gunzip -c \$genome_fp | sed 's/ .*//' | gzip -c >> ALL_GENOMES.fasta.gz

    # Make a headerless CSV with the first column containing the
    # headers from this FASTA, and the second containing a name
    # for this organism which includes both the #Organism Name
    # as well as the Assembly

    # Concatenate this text to the file ALL_GENOMES.csv.gz

    gunzip -c \$genome_fp | grep '>' | tr -d '>' | while read header; do
        echo \$header,\$genome_name
    done | gzip -c >> ALL_GENOMES.csv.gz
    

    echo "\\n\\n"
done

    """

}

workflow {

    // Read four columns from the manifest CSV
    //  - #Organism Name
    //  - Assembly
    //  - GenBank FTP

    // Skip any rows that are lacking the GenBank FTP column

    Channel.from(
        file(params.manifest)
    ).splitCsv(
        header: true
    ).filter {
        r -> r["GenBank FTP"].length() > 0
    }.map {
        r -> "${r['#Organism Name'].replaceAll(/|/, '').replaceAll(/ /, '_')} ${r['Assembly'].replaceAll(/ /, '_')} ${r['GenBank FTP'].replaceAll(/"/, "").replaceAll(/ /, '_')}"
    }.toSortedList(
    ).flatten(
    ).collate(
        params.batchsize
    ).map { 
        it -> it.join("\n") 
    }.set {
        genome_ch
    }

    // Download the genomes
    downloadGenome(
        genome_ch
    )

    // Concatenate the genomes in two steps
    concatGenomes_1(
        downloadGenome.out[0].toSortedList().flatten().collate(100),
        "concat_genomes_1.fasta.gz"
    )
    concatGenomes_2(
        concatGenomes_1.out.collect(),
        "${params.dbname}.fasta.gz"
    )
    // Concatenate the CSV header key in two steps
    concatCSV_1(
        downloadGenome.out[1].toSortedList().flatten().collate(100),
        "concat_genomes_1.csv.gz"
    )
    concatCSV_2(
        concatCSV_1.out.collect(),
        "${params.dbname}.csv.gz"
    )

    publish:

    concatGenomes_2.out to: "${output_folder}", mode: "copy"
    concatCSV_2.out to: "${output_folder}", mode: "copy"
}