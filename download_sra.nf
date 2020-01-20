#!/usr/bin/env nextflow

// Script to download FASTQ files from SRA

// Set default values for parameters
params.accession = false
params.output = false
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/geneshot/download_sra <ARGUMENTS>

    NOTE:   This script expected paired-end FASTQ data, and will probably not work
            with any other type
    
    Required Arguments:
      --accession           Accession for SRA entries to download, comma separated
      --output              Folder to write output files to

    Output Files:

    All output files will be written to the --output folder. This includes one or two
    FASTQ files per Run as well as a manifest.csv file listing all of the files which
    were downloaded.

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || params.accession == false || params.output == false){
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

// Make a channel with the files needed for every Run
process downloadSRA {
    container "quay.io/fhcrc-microbiome/get_sra@sha256:16b7988e435da5d21bb1fbd7c83e97db769f1c95c9d32823fde49c729a64774a"
    label "io_limited"
    errorStrategy 'retry'
    publishDir "${output_folder}", mode: "copy"
    
    input:
    val accession from Channel.from(params.accession.split(","))

    output:
    tuple accession, file("*fastq.gz") into reads_ch


    """
    set -e

    echo "Setting up the cache folder"
    mkdir cache
    vdb-config --root -s /repository/user/main/public/root=\$PWD/cache || echo Setting up the download cache

    accession=$accession
    echo "Downloading \$accession"

    fastq-dump \
        --split-files \
        --outdir ./ \
        --readids \
        \$accession

    echo "Compressing downloaded FASTQ files"
    gzip \$accession*

    echo "Removing the cache"
    rm -rf cache

    echo "Done"
    """

}

// Make a comma-separated string with all of the files which were downloaded
manifestStr = reads_ch.reduce(
    'specimen,R1,R2\n'
){ csvStr, row ->
    return  csvStr += "${row[0]},${output_folder}${row[1][0].name},${output_folder}${row[1][1].name}\n";
}


process outputManifest {
    container "ubuntu:18.04"

    publishDir path: "${output_folder}", mode: "copy"

    input:
        val manifestStr
    
    output:
        file 'manifest.csv'

"""
#!/bin/bash

echo "${manifestStr}" > manifest.csv
"""
}
