#!/usr/bin/env nextflow

// Default values for boolean flags
params.interleaved = false
params.paired = false
params.help = false

def helpMessage() {
    log.info"""
    Usage:

    nextflow run fredhutch/geneshot <ARGUMENTS>
    
    Arguments:
      --batchfile                   CSV file listing samples to analyze (see below)
      --ref_dmnd                    Path to reference database in DIAMOND (.dmnd) format
      --ref_hdf5                    Path to HDF5 file containing reference database metadata
      --output_folder               Folder to place outputs
      --output_prefix               Name for output files

    Options:
      --paired                      Input data is paired-end FASTQ in two files (otherwise treat as single-ended)
      --interleaved                 Input data is interleaved paired-end FASTQ in one file (otherwise treat as single-ended)

    Batchfile:
      The batchfile is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `name`. 
      Default is to expect a single column `fastq` pointing to a single-ended or interleaved FASTQ.
      If data is --paired, reads are specified by two columns, `fastq1` and `fastq2`.

    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
if (params.help){
    helpMessage()
    exit 0
}
// --output_folder is the folder in which to place the results
params.output_folder = "./"

// --output_prefix is the name to prepend to all output files
params.output_prefix = "geneshot_output"

// Logic to handle different types of input data
if ( params.paired ){
  assert !params.interleaved: "--paired cannot be specified together with --interleaved"

  Channel.from(file(params.batchfile))
        .splitCsv(header: true, sep: ",")
        .map { sample ->
        [sample.name, file(sample.fastq1), file(sample.fastq2)]}
        .set{ interleave_ch }

  process interleave {
    container "ubuntu:16.04"
    cpus 4
    memory "8 GB"
    errorStrategy "retry"

    input:
    set sample_name, file(fastq1), file(fastq2) from interleave_ch

    output:
    set sample_name, file("${fastq1}.interleaved.fastq.gz") into concatenate_ch

    afterScript "rm *"

    """
    set -e

    # Some basic checks that the files exist and the line numbers match
    [[ -s "${fastq1}" ]]
    [[ -s "${fastq2}" ]]
    (( \$(gunzip -c ${fastq1} | wc -l) == \$(gunzip -c ${fastq2} | wc -l) ))

    # Now interleave the files
    paste <(gunzip -c ${fastq1}) <(gunzip -c ${fastq2}) | paste - - - - | awk -v OFS="\\n" -v FS="\\t" '{print(\$1,\$3,\$5,\$7,\$2,\$4,\$6,\$8)}' | gzip -c > "${fastq1}.interleaved.fastq.gz"
    """
      
  }

}
else {
    // For either `interleaved` or `single` (no flags), just put the FASTQ into the same analysis queue
    Channel.from(file(params.batchfile))
        .splitCsv(header: true, sep: ",")
        .map { sample ->
        [sample.name, file(sample.fastq)]}
        .set{ concatenate_ch }

}

// Concatenate reads by sample name
process concatenate {
  container "ubuntu:16.04"
  cpus 4
  memory "8 GB"
  errorStrategy "retry"
  
  input:
  set sample_name, file(fastq_list) from concatenate_ch.groupTuple()
  
  output:
  set sample_name, file("${sample_name}.fastq.gz") into count_reads, metaphlan_ch, diamond_ch

  afterScript "rm *"

  """
set -e
ls -lahtr
cat ${fastq_list} > TEMP && mv TEMP ${sample_name}.fastq.gz
  """

}

// Count the number of input reads
process countReads {
  container "ubuntu:16.04"
  cpus 1
  memory "4 GB"
  
  input:
  set sample_name, file(fastq) from count_reads
  
  output:
  file "${sample_name}.countReads.csv" into total_counts

  afterScript "rm *"

  """
set -e
n=\$(gunzip -c "${fastq}" | awk 'NR % 4 == 1' | wc -l)
echo "${sample_name},total_reads,\$n" > "${sample_name}.countReads.csv"
  """

}

process metaphlan2 {
    container "quay.io/fhcrc-microbiome/metaphlan@sha256:51b416458088e83d0bd8d840a5a74fb75066b2435d189c5e9036277d2409d7ea"
    cpus 16
    memory "32 GB"

    input:
    set val(sample_name), file(input_fastq) from metaphlan_ch
    
    output:
    file "${sample_name}.metaphlan.tsv" into metaphlan_for_summary

    afterScript "rm *"

    """
    set -e
    metaphlan2.py --input_type fastq --tmp_dir ./ -o ${sample_name}.metaphlan.tsv ${input_fastq}
    """
}

process diamond {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    cpus 32
    memory "240 GB"
    errorStrategy "retry"
    
    input:
    set val(sample_name), file(input_fastq) from diamond_ch
    file refdb from file(params.ref_dmnd)
    val min_id from 90
    val query_cover from 50
    val cpu from 32
    val top from 1
    val min_score from 20
    val blocks from 15
    val query_gencode from 11

    output:
    set sample_name, file("${sample_name}.aln.gz") into aln_ch

    afterScript "rm *"

    """
    set -e
    diamond \
      blastx \
      --query ${input_fastq} \
      --out ${sample_name}.aln.gz \
      --threads ${cpu} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --min-score ${min_score} \
      --query-cover ${query_cover} \
      --id ${min_id} \
      --top ${top} \
      --block-size ${blocks} \
      --query-gencode ${query_gencode} \
      --compress 1 \
      --unal 0
    """

}


process famli {
    container "quay.io/fhcrc-microbiome/famli@sha256:241a7db60cb735abd59f4829e8ddda0451622b6eb2321f176fd9d76297d8c9e7"
    cpus 16
    memory "120 GB"
    errorStrategy "retry"
    
    input:
    set sample_name, file(input_aln) from aln_ch
    val cpu from 16
    val batchsize from 50000000

    output:
    file "${sample_name}.json.gz" into famli_json_for_summary

    afterScript "rm *"

    """
    set -e
    famli \
      filter \
      --input ${input_aln} \
      --output ${sample_name}.json \
      --threads ${cpu} \
      --batchsize ${batchsize}
    gzip ${sample_name}.json
    """

}

process summarizeExperiment {
    container "quay.io/fhcrc-microbiome/python-pandas:latest"
    cpus 16
    memory "32 GB"
    publishDir "${params.output_folder}"

    input:
    file metaphlan_tsv_list from metaphlan_for_summary.collect()
    file famli_json_list from famli_json_for_summary.collect()
    file ref_hdf5 from file(params.ref_hdf5)
    val output_prefix from params.output_prefix

    output:
    file "${output_prefix}.hdf5"

    afterScript "rm *"

"""
#!/usr/bin/env python3
import gzip
import json
import logging
import os
import pandas as pd

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [assembleAbundances] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write logs to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)

# Read in all of the FAMLI results
allele_abund = {}
for fp in os.listdir("./"):
    if fp.endswith(".json.gz"):
        sample_name = fp.replace(".json.gz", "")
        print("Reading in %s" % (fp))
        allele_abund[sample_name] = [
            dict([
                ("id": r["id"]),
                ("depth": r["depth"])
            ])
            for r in json.load(gzip.open(fp, "rt"))
        ]

"""

}