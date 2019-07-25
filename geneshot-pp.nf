#!/usr/bin/env nextflow

// Geneshot preprocessing submodule

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.index = false
params.help = false
params.adapter_F = "CTGTCTCTTATACACATCT"
params.adapter_R = "CTGTCTCTTATACACATCT"

params.output_folder = '.'


// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run geneshot-pp.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest         CSV file listing samples to preprocess (see below)
      --output_folder    Folder to place outputs
      --output_prefix    Name for output files

    Options:
      --index            Index reads are provided (default: false)

    Batchfile:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `name`. 
      Reads are specified by two columns, `fastq1` and `fastq2`.
      If index reads are provided, the column titles should be 'index1' and 'index2'

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

if (params.index) {
    Channel.from(file(params.manifest))
        .splitCsv(header: true, sep: ",")
        .map { sample ->
        [sample.name, file(sample.fastq1), file(sample.fastq2), file(sample.index1), file(sample.index2)]}
        .set{ input_ch }
    // implement barcodecop here
}
else {
    Channel.from(file(params.manifest))
        .splitCsv(header: true, sep: ",")
        .map { sample ->
        [sample.name, file(sample.fastq1), file(sample.fastq2)]}
        .set{ demupltiplexed_ch }
}

process cutadapt {
  container "golob/cutadapt:1.18__bcw.0.3.0_al38"
  cpus 1
  memory "4 GB"
  errorStrategy "retry"

  publishDir "${params.output_folder}/noadapt/"

  input:
  set sample_name, file(fastq1), file(fastq2) from demupltiplexed_ch
  
  output:
  set sample_name, file("${fastq1}.noadapt.fastq.gz"), file("${fastq2}.noadapt.fastq.gz"), file("${fastq1}.cutadapt.log") into noadapt_ch

  """
  cutadapt -a \
  ${params.adapter_F} -A ${params.adapter_R} \
  -o ${fastq1}.noadapt.fastq.gz -p ${fastq2}.noadapt.fastq.gz \
  ${fastq1} ${fastq2} > ${fastq1}.cutadapt.log
  """
}