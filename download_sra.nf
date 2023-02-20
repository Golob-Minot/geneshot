#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Script to download FASTQ files from SRA

// Set default values for parameters
params.accession = false
params.list = false
params.output = '.'
params.help = false

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/geneshot/download_sra.nf <ARGUMENTS>

    NOTE:   This script expects paired-end FASTQ data, and will not download any other type

    Required Arguments:
      --accession           Accession for NCBI BioProject to download
                    OR
      --list                List of SRA IDs to download
      --output              Folder to write output files (default invoking directory)


    Output Files:

    All output files will be written to the --output folder. This includes one or two
    FASTQ files per Run as well as a `BIOPROJECT.metadata.csv` file listing all of the files which
    were downloaded, as well as the metadata describing those samples within NCBI.

    The `BIOPROJECT.metadata.csv` file will also include the metadata recorded for this set of Runs
    within the SRA database. The columns for this file may not be formatted nicely,
    but they do match the structure of the data within the SRA API.

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
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

// Main workflow
workflow {

    if (params.accession) {
        // Get a file listing every SRA ID
        getSRAlist(
            params.accession
        )

        // Make a channel with every individual SRA ID
        sra_id_ch = getSRAlist.out.splitText().map { r -> r.replaceAll(/\n/, "")}
    } else if (params.list) {
        sra_acc_ch = Channel.fromList(
            file(params.list).readLines().each { it.strip() }
       )

    
        sra_id_ch = getSRRuid(
            sra_acc_ch 
        ).map { it.strip() }
        
    } else {
        helpMessage()
        exit 0
    }

    // Get the metadata for each accession
    getMetadata(
        sra_id_ch
    )

    // Join together all of the metadata into a single table
//    joinMetadata(
//        getMetadata.out.toSortedList()
//    )

    // Parse the accession from the metadata file names
    sra_acc_ch = getMetadata.out.map {
        r -> r.name.replaceAll(/.metadata.json.gz/, "")
    }

    // Fetch the SRA files for each accession
    downloadSRA(
        sra_acc_ch
    )

    // Extract the FASTQ files for each accession
    extractSRA(
        downloadSRA.out
    )

    // Make a comma-separated string with all of the files which were downloaded
    manifestStr = extractSRA.out.reduce(
        'specimen,R1,R2\n'
    ){ csvStr, row ->
        return  csvStr += "${row[0]},${output_folder}${row[1][0].name},${output_folder}${row[1][1].name}\n";
    }

    // Make a single set of readnames joined with the metadata
//    gatherReadnames(
//        manifestStr,
//        joinMetadata.out
//    )

}

// Get the accession for each Run in this BioProject
process getSRAlist {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_net"
    errorStrategy 'finish'
    
    input:
    val accession

    output:
    file "accession_list.txt"

"""
#!/usr/bin/env python3
import requests
from time import sleep
import xml.etree.ElementTree as ET
import pandas as pd

accession ="${accession}"
print("Fetching paired-end data from %s" % accession)

# Make a function to fetch data from the Entrez API
# while retrying if any errors are encountered
def request_url(url, max_retries=10, pause_seconds=0.5):
    # Try the request `max_retries` times, with a `pause_seconds` pause in-between
    assert isinstance(max_retries, int)
    assert max_retries > 0
    for i in range(max_retries):
        r = requests.get(url)
        if r.status_code == 200:
            return r.text
        print("Caught an error on attempt #%d, retrying" % (i + 1))
        sleep(pause_seconds)

# Format an Entrez query, execute the request, and return the result
def entrez(mode, **kwargs):
    assert mode in ["efetch", "esearch", "elink"]

    kwargs_str = "&".join(["%s=%s" % (k, v) for k, v in kwargs.items()])
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/%s.fcgi?%s"
    url = url % (mode, kwargs_str)

    # Get the text response
    response_text = request_url(url)

    # Parse the XML
    return ET.fromstring(response_text)

# Make a list of SRA records which will be added to
sra_id_list = []

# First, get the ID(s) for the BioProject from the accession
bioproject_result = entrez('esearch', db="bioproject", term=accession)

# Iterate over each of the IDs linked to the BioProject
for i in bioproject_result.find(
    "IdList"
).findall(
    "Id"
):
    # Looking at this particular ID
    bioproject_id = i.text

    print("Getting BioProject %s = %s" % (accession, bioproject_id))

    # Get all of the SRA Runs for this BioProject
    elink_results = entrez(
        'elink',
        dbfrom="bioproject", 
        id=int(bioproject_id), 
        linkname="bioproject_sra"
    )

    # Make sure that there are some links in this set of results, or move on
    if elink_results.find("LinkSet") is None:
        print("No SRA accessions found for %s = %s" % (accession, bioproject_id))
        continue
    if elink_results.find("LinkSet").find("LinkSetDb") is None:
        print("No SRA accessions found for %s = %s" % (accession, bioproject_id))
        continue

    # Parse all of the SRA records from this BioProject
    sra_id_list.extend([
        child.find("Id").text
        for child in elink_results.find("LinkSet").find("LinkSetDb")
        if child.find("Id") is not None
    ])

print("Found %d SRA accessions from this BioProject" % (len(sra_id_list)))

assert len(sra_id_list) > 0

# Save the accession list
with open("accession_list.txt", "wt") as fo:
    fo.write("\\n".join(sra_id_list))

"""
}

// Get the UID for an SRR
process getSRRuid {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_net"
    errorStrategy 'ignore'
    
    input:
    val accession

    output:
    stdout

"""
#!/usr/bin/env python3
import requests
from time import sleep
import xml.etree.ElementTree as ET
import os

run_accession ="${accession}"

# Make a function to fetch data from the Entrez API
# while retrying if any errors are encountered
def request_url(url, max_retries=10, pause_seconds=0.5):
    # Try the request `max_retries` times, with a `pause_seconds` pause in-between
    assert isinstance(max_retries, int)
    assert max_retries > 0
    for i in range(max_retries):
        r = requests.get(url)
        if r.status_code == 200:
            return r.text
        sleep(pause_seconds)

# Format an Entrez query, execute the request, and return the result
def entrez(mode, **kwargs):
    assert mode in ["efetch", "esearch", "elink"]

    kwargs_str = "&".join(["%s=%s" % (k, v) for k, v in kwargs.items()])
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/%s.fcgi?%s"
    url = url % (mode, kwargs_str)

    # Get the text response
    response_text = request_url(url)

    # Parse the XML
    return ET.fromstring(response_text)


# First, serch for the SRA accession
SRA_result = entrez(
    'esearch', 
    db="sra", 
    term=run_accession
)

SRR_id = SRA_result.find(
    "IdList"
).find(
    "Id"
).text

print(SRR_id.strip())

"""
}


// Get the metadata for a single SRA accession
process getMetadata {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_net"
    errorStrategy 'ignore'
    
    input:
    val sra_id

    output:
    file "*.metadata.json.gz"

"""
#!/usr/bin/env python3
import json
import gzip
import requests
from time import sleep
import xml.etree.ElementTree as ET
import pandas as pd

print("Fetching metadata for SRA ID ${sra_id}")

# Make a function to fetch data from the Entrez API
# while retrying if any errors are encountered
def request_url(url, max_retries=10, pause_seconds=0.5):
    # Try the request `max_retries` times, with a `pause_seconds` pause in-between
    assert isinstance(max_retries, int)
    assert max_retries > 0
    for i in range(max_retries):
        r = requests.get(url)
        if r.status_code == 200:
            return r.text
        print("Caught an error on attempt #%d, retrying" % (i + 1))
        sleep(pause_seconds)

# Format an Entrez query, execute the request, and return the result
def entrez(mode, **kwargs):
    assert mode in ["efetch", "esearch", "elink"]

    kwargs_str = "&".join(["%s=%s" % (k, v) for k, v in kwargs.items()])
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/%s.fcgi?%s"
    url = url % (mode, kwargs_str)

    # Get the text response
    response_text = request_url(url)

    # Parse the XML
    return ET.fromstring(response_text)

# Function to iteratively parse the XML record
def parse_record(r, prefix=[]):

    if r.text is not None:
        if len(r.text.replace("\\n", "")) > 0:
            yield prefix + [r.tag], r.text.replace("\\n", "")

    for k, v in r.attrib.items():
        yield prefix + [k], v

    for child in r:

        if child.tag == "TAG":
            if r.find("VALUE") is not None:
                yield prefix + [child.text], r.find("VALUE").text
                continue

        if child.tag == "VALUE":
            if r.find("TAG") is not None:
                continue

        for i in parse_record(child, prefix + [r.tag]):
            yield i

# Function to get metadata from an SRA accession
def fetch_sra_metadata(sra_id):
    print("Fetching metadata for %s" % sra_id)
    record = entrez(
        'efetch',
        db="sra",
        id=sra_id
    )

    # Keep track of the data as a simple dict
    dat = {}

    dat["LAYOUT"] = record.find("EXPERIMENT_PACKAGE").find("EXPERIMENT").find("DESIGN").find("LIBRARY_DESCRIPTOR").find("LIBRARY_LAYOUT")[0].tag

    # Iterate over every entry in the record
    for prefix, val in parse_record(record):

        # Lots of things to skip as not being worth saving
        if len(prefix) < 2:
            continue
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


    # Make sure that we have the RUN_SET_accession
    assert "RUN_SET_accession" in dat

    return dat

# Fetch the metadata
dat = fetch_sra_metadata("${sra_id}")

# Get the SRR accession from the ID that we used as the input
sra_accession = dat["RUN_SET_accession"]

# Write out to a file
with gzip.open("%s.metadata.json.gz" % sra_accession, "wt") as fo:
    fo.write(
        json.dumps(
            dat,
            indent=4
        )
    )

"""
}

// Collect the metadata from all accessions
process joinMetadata {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "mem_medium"
    errorStrategy 'ignore'
    
    input:
    file metadata_json_list

    output:
    file "${params.accession}.metadata.csv"

"""
#!/usr/bin/env python3
import os
import json
import gzip
import pandas as pd

# Parse the list of files to read in
fp_list = []
for fp in os.listdir("."):
    if ".metadata.json.gz" in fp:
        assert fp.endswith(".metadata.json.gz")
        fp_list.append(fp)
print("Reading in a set of %d metadata files" % len(fp_list))

# Iterate over all of the SRA IDs from this BioProject
metadata_df = pd.DataFrame([
    json.load(gzip.open(fp, "rt"))
    for fp in fp_list
])
print("Found records for %d SRA accessions" % metadata_df.shape[0])

# Make sure that we have a column for the SRR* accession
assert "RUN_SET_accession" in metadata_df.columns.values
assert metadata_df["RUN_SET_accession"].isnull().sum() == 0

# Save the complete metadata to a file
metadata_df.to_csv("${params.accession}.metadata.csv", index=None)

"""
}

// Download the .sra file for each SRR accession
process downloadSRA {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_net"
    errorStrategy 'ignore'
    
    input:
    val accession

    output:
    path "${accession}", optional: true


"""
set -e

echo "Getting the URL for the SRA file"
ACC=$accession
curl -o \${ACC}.json -s -X POST "https://www.ncbi.nlm.nih.gov/Traces/sdl/2/retrieve?acc=\${ACC}"

sra_url="\$(cat \${ACC}.json | jq '.result | .[0] | .files | .[0] | .locations | .[0] | .link' | tr -d '"')"
echo "Download URL is \$sra_url"

if [[ \$sra_url == null ]]; then
    cat \${ACC}.json
    echo "Stopping"
else
    echo "Downloading"
    wget -O \${ACC} \${sra_url}

    echo "Done"
fi
"""

}

// Extract the FASTQ files from the SRA file
process extractSRA {
    container "golob/get_sra:v3.0.0A"
    label "io_net"
    errorStrategy 'ignore'
    publishDir "${output_folder}", mode: "copy", overwrite: "true"
    
    input:
    file accession

    output:
    tuple val("${accession.name}"), file("*fastq.gz")


    """
    set -e

    fastq-dump \
        --split-files \
        --outdir ./ \
        ${accession}

    rm ${accession}

    echo "Compressing downloaded FASTQ files"
    pigz ${accession.name}*

    echo "Done"
    """

}


process gatherReadnames {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label "io_limited"
    errorStrategy 'finish'
    publishDir "${output_folder}", mode: "copy", overwrite: "true"

    input:
        val manifestStr
        file metadata_csv
    
    output:
        file "${params.accession}.csv"

"""
#!/usr/bin/env python3
import pandas as pd

with open("manifest.csv", "wt") as fo:
    fo.write(\"\"\"${manifestStr}\"\"\")

manifest_df = pd.read_csv("manifest.csv")
print(manifest_df)

# Read in the metadata
metadata_df = pd.read_csv("${metadata_csv}")
print(metadata_df)

# Join the two together
metadata_df = pd.concat([
    manifest_df.set_index("specimen"),
    metadata_df.set_index("RUN_SET_accession")
], axis=1, sort=True)

# Write out to a file
metadata_df.reset_index(
).rename(
    columns=dict([("index", "specimen")])
).to_csv(
    "${params.accession}.csv", 
    index=None
)
"""
}
