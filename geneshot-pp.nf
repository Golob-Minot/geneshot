#!/usr/bin/env nextflow

/*
  Geneshot preprocessing submodule:
    Steps:
    1) (if index is available): barcodecop to verify demultiplexing
    2) cutadapt to remove adapters.
    3) remove human reads via
      3A) downloading the cached human genome index
      3B) aligning against the human genome and extracting unpaired reads
*/

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.index = false
params.help = false
params.adapter_F = "CTGTCTCTTATACACATCT"
params.adapter_R = "CTGTCTCTTATACACATCT"
params.hg_index_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30

params.output_folder = '.'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run geneshot-pp.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples to preprocess (see below)

    Options:
      --output_folder       Folder to place outputs (default invocation dir)
      --index               Index reads are provided (default: false)
      --hg_index_url        URL for human genome index, defaults to current HG
      --hg_index            Cached copy of the bwa indexed human genome, TGZ format
      --adapter_F           Forward sequencing adapter sequence (to be removed)
      --adapter_R           Reverse sequencing adapter sequence (to be removed)
                              (Adapter sequences default to nextera adapters)
      --min_hg_align_score  Minimum alignment score for human genome (default 30)
      -w                    Working directory. Defaults to `./work`

    Batchfile:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This can be repeated. 
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If index reads are provided, the column titles should be 'I1' and 'I2'

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
        [sample.specimen, file(sample.R1), file(sample.R2), file(sample.I1), file(sample.I2)]}
        .set{ input_ch }
    // Step 1: barcodecop
    process barcodecop {
      container "golob/barcodecop:0.4.1__bcw_0.3.0"
      label 'io_limited'
      errorStrategy "retry"

      input:
        set specimen, file(R1), file(R2), file(I1), file(I2) from input_ch
      
      output:
        set specimen, file("${R1}.bcc.fq.gz"), file("${R2}.bcc.fq.gz") into demupltiplexed_ch
      """
      set -e

      barcodecop \
      ${I1} ${I2} \
      --match-filter \
      -f ${R1} \
      -o ${R1}.bcc.fq.gz &&
      barcodecop \
      ${I1} ${I2} \
      --match-filter \
      -f ${R2} \
      -o ${R2}.bcc.fq.gz
      """
    }
}
else {
    Channel.from(file(params.manifest))
        .splitCsv(header: true, sep: ",")
        .map { sample ->
        [sample.specimen, file(sample.R1), file(sample.R2)]}
        .set{ demupltiplexed_ch }
}

// Step 2
process cutadapt {
  container "golob/cutadapt:2.3__bcw.0.3.0_al38B_FH"
  label 'io_limited'
  errorStrategy "retry"

  //publishDir "${params.output_folder}/noadapt/"

  input:
  set sample_name, file(fastq1), file(fastq2) from demupltiplexed_ch
  
  output:
  set sample_name, file("${fastq1}.noadapt.R1.fq.gz"), file("${fastq2}.noadapt.R2.fq.gz"), file("${fastq1}.cutadapt.log") into noadapt_ch

  """
  set -e 

  cutadapt \
  -j ${task.cpus} \
   -a ${params.adapter_F} -A ${params.adapter_R} \
  -o ${fastq1}.noadapt.R1.fq.gz -p ${fastq2}.noadapt.R2.fq.gz \
  ${fastq1} ${fastq2} > ${fastq1}.cutadapt.log
  """
}

// Step 3A.
if (!params.hg_index) {
  process download_hg_index {
    container "golob/bwa:0.7.17__bcw.0.3.0I"
    errorStrategy "retry"
    label 'io_limited'

    output:
      file 'hg_index.tar.gz' into hg_index_tgz
    
    """
    wget ${params.hg_index_url} -O hg_index.tar.gz
    """
  }
} else {
  hg_index_tgz = file(params.hg_index)
}


// Step 3B.
process remove_human {
  container "golob/bwa:0.7.17__bcw.0.3.0I"
  //container "quay.io/fhcrc-microbiome/bwa:v0.7.17--4"
  errorStrategy "retry"
  publishDir "${params.output_folder}/nohuman/"
  label 'mem_veryhigh'


  input:
    file hg_index_tgz from hg_index_tgz
    set sample_name, file(fastq1), file(fastq2), file(cutadapt_log) from noadapt_ch
  
  output:
    set sample_name, file("${fastq1}.nohuman.R1.fq.gz"), file("${fastq2}.nohuman.R2.fq.gz"), file("${fastq1}.nohuman.log") into nohuman_ch

  afterScript "rm -rf hg_index/*"

  """
  set - e

  bwa_index_prefix=\$(tar -ztvf ${hg_index_tgz} | head -1 | sed \'s/.* //\' | sed \'s/.amb//\') && \
  echo BWA index file prefix is \${bwa_index_prefix} | tee -a ${fastq1}.nohuman.log && \
  echo Extracting BWA index | tee -a ${fastq1}.nohuman.log && \
  mkdir -p hg_index/ && \
  tar -I pigz -xf ${hg_index_tgz} -C hg_index/ | tee -a ${fastq1}.nohuman.log && \
  echo Files in index directory: | tee -a ${fastq1}.nohuman.log && \
  ls -l -h hg_index | tee -a ${fastq1}.nohuman.log && \
  echo Running BWA | tee -a ${fastq1}.nohuman.log && \
  bwa mem -t ${task.cpus} \
  -T ${params.min_hg_align_score} \
  -o alignment.sam \
  hg_index/\$bwa_index_prefix \
  ${fastq1} ${fastq2} \
  | tee -a ${fastq1}.nohuman.log && \
  echo Checking if alignment is empty  | tee -a ${fastq1}.nohuman.log && \
  [[ -s alignment.sam ]] && \
  echo Extracting Unaligned Pairs | tee -a ${fastq1}.nohuman.log && \
  samtools fastq alignment.sam \
  --threads ${task.cpus} -f 12 \
  -1 ${fastq1}.nohuman.R1.fq.gz -2 ${fastq2}.nohuman.R2.fq.gz \
  | tee -a ${fastq1}.nohuman.log && \
  echo Done | tee -a ${fastq1}.nohuman.log
  """
}

//   gunzip -c ${hg_index_tgz} | tar -x -C hg_index/ | tee -a ${fastq1}.nohuman.log && \
//  tar -I pigz -xf ${hg_index_tgz} -C hg_index/ | tee -a ${fastq1}.nohuman.log && \

nohuman_ch.reduce('specimen, R1, R2\n'){ csvStr, row ->
            return  csvStr += "${row[0]}, ${row[1]}, ${row[2]}\n";
        }.set{manifestStr}

process outputManifest {
    container "ubuntu:16.04"

    publishDir "${params.output_folder}/"

    input:
        val manifestStr from manifestStr
    
    output:
        file 'manifest.nohuman.csv' into opManifest

    """
        echo "${manifestStr}" > manifest.nohuman.csv
    """
}
