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
      --accession           Accession for NCBI BioProject to download
      --output              Folder to write output files

    Output Files:

    All output files will be written to the --output folder. This includes one or two
    FASTQ files per Run as well as a manifest.csv file listing all of the files which
    were downloaded.

    The manifest.csv file will also include the metadata recorded for this set of Runs
    within the SRA database. The columns for this file may not be formatted nicely,
    but they do match the structure of the data within the SRA API.

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

// Get the accession for each Run in this BioProject
process getRunAccessions {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_limited"
    errorStrategy 'retry'
    
    input:
    val accession from params.accession

    output:
    file "accession_list.txt" into accession_list


"""
#!/usr/bin/env python3
from Bio import Entrez
from time import sleep
import xml.etree.ElementTree as ET

Entrez.email = "scicomp@fredhutch.org"

accession ="${accession}"
print(accession)

# First, get the ID for the BioProject from the accession
try:
    handle = Entrez.esearch(db="bioproject", term=accession)
except:
    sleep(0.2)
    handle = Entrez.esearch(db="bioproject", term=accession)
l = "".join(handle.readlines())
handle.close()

bioproject_id = None
for i in l.split("<Id>"):
    if "</Id>" in i:
        bioproject_id = i.split("</Id>")[0]
assert bioproject_id is not None

print("Getting BioSample %s = %s" % (accession, bioproject_id))

# Get all of the SRA Runs for this BioProject
try:
    handle = Entrez.elink(
        dbfrom="bioproject", 
        id=int(bioproject_id), 
        linkname="bioproject_sra"
    )
except:
    sleep(0.2)
    handle = Entrez.elink(
        dbfrom="bioproject", 
        id=int(bioproject_id), 
        linkname="bioproject_sra"
    )
l = "".join(handle.readlines())
handle.close()

with open("accession_list.txt", "wt") as handle:
    for i in "".join(l).split("<Id>"):
        if "</Id>" in i:
            if bioproject_id in i:
                continue
            handle.write(
                i.split("</Id>")[0] + "\\n"
            )
    
"""
}

// Get the metadata for each SRA accession, including the name used to download the FASTQ
process getMetadata {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_limited"
    errorStrategy 'retry'
    
    input:
    val accession from accession_list.splitText()

    output:
    file "*.csv" into metadata_csv_ch1
    file "*.csv" into metadata_csv_ch2


"""
#!/usr/bin/env python3
from Bio import Entrez
import xml.etree.ElementTree as ET
import pandas as pd
from time import sleep
Entrez.email = "scicomp@fredhutch.org"

accession = "${accession.replaceAll(/\n/, "")}"
print(accession)

# Function to iteratively parse the XML record
def parse_record(r, prefix=[]):
    for k, v in r.attrib.items():
        yield prefix + [k], v
    for child in r:
        for i in parse_record(child, prefix + [r.tag]):
            yield i

try:
    handle = Entrez.efetch(db="sra", id=int(accession))
except:
    sleep(1)
    handle = Entrez.efetch(db="sra", id=int(accession))

record = ET.fromstring("".join(handle.readlines()))
handle.close()

# Keep track of the data as a simple dict
dat = {}

# Iterate over every entry in the record
for prefix, val in parse_record(record):

    # Lots of things to skip as not being worth saving
    if prefix[0] != "EXPERIMENT_PACKAGE_SET":
        continue
    if prefix[1] != "EXPERIMENT_PACKAGE":
        continue
    if len(val) == 0:
        continue
    if prefix[-1] == "namespace":
        continue
    if "SRAFiles" in prefix:
        continue
    if "CloudFiles" in prefix:
        continue
    if "Statistics" in prefix:
        continue
    if "Databases" in prefix:
        continue

    # This item is worth saving
    dat["_".join(prefix[2:])] = val

# Get the accession from this set of data
assert "RUN_SET_accession" in dat
accession = dat["RUN_SET_accession"]

# Write out the data as a simple, one-line CSV
pd.DataFrame([dat]).to_csv(accession + ".csv", index=None)

"""
}

// Make a channel with the files needed for every Run
process downloadSRA {
    container "quay.io/fhcrc-microbiome/get_sra@sha256:16b7988e435da5d21bb1fbd7c83e97db769f1c95c9d32823fde49c729a64774a"
    label "io_limited"
    errorStrategy 'retry'
    publishDir "${output_folder}", mode: "copy"
    
    input:
    val metadata_csv from metadata_csv_ch1

    output:
    tuple val("${metadata_csv.name.replaceAll(/.csv/, "")}"), file("*fastq.gz") into reads_ch


    """
    set -e

    echo "Setting up the cache folder"
    mkdir cache
    vdb-config --root -s /repository/user/main/public/root=\$PWD/cache || echo Setting up the download cache

    accession=${metadata_csv.name.replaceAll(/.csv/, "")}
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


process gatherReadnames {
    container "ubuntu:18.04"

    input:
        val manifestStr
    
    output:
        file 'manifest.csv' into manifest_csv

"""
#!/bin/bash

echo "${manifestStr}" > manifest.csv
"""
}


process gatherMetadata {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_limited"
    errorStrategy 'retry'
    publishDir "${output_folder}", mode: "copy"
    
    input:
    file csv_list from metadata_csv_ch2.toSortedList()
    file "input_manifest.csv" from manifest_csv

    output:
    file "manifest.csv"


"""
#!/usr/bin/env python3
import pandas as pd
import os

# Read in the manifest with the location of R1 and R2
manifest = pd.read_csv(
    "input_manifest.csv"
).set_index(
    "specimen"
)

csv_list = "${csv_list}".split(" ")

# Read in all of the CSV files
df = {}

# Iterate over every file in the list
for fp in csv_list:

    # Make sure that the file name ends with ".csv" as expected
    assert fp.endswith(".csv"), fp

    # Get the accession from the file name
    accession = fp.replace(".csv", "")

    # Read in the CSV
    df[accession] = pd.read_csv(fp).iloc[0]

# Make a single DataFrame
df = pd.DataFrame(
    df
).T

# Join together the two tables
for k in df.columns.values:
    if k not in manifest.columns.values:
        print(k)
        manifest[k] = df[k]

# Write out to the file
manifest.reset_index().to_csv("manifest.csv", index=None)
"""

}
