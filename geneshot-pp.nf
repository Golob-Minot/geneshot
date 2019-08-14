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
if (params.help || params.manifest == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

// Implicit else we have a manifest. Time to do some validation:
// Split into index provided or not. Further validate critical elements are present
// e.g. A specimen, R1 and R2

input_w_index_ch = Channel.create()
input_no_index_ch = Channel.create()
input_invalid_ch = Channel.create()

Channel.from(file(params.manifest))
    .splitCsv(header: true, sep: ",")
    .choice(
        input_invalid_ch,
        input_no_index_ch,
        input_w_index_ch,
    ) {
        r -> if (
            (r.specimen == null) ||
            (r.R1 == null) ||
            (r.R2 == null) ||
            (r.specimen == "") ||
            (r.R1 == "") ||
            (r.R2 == "")
        ) return 0;
        else if (
            (r.I1 != null) && 
            (r.I2 != null) && 
            (r.I1 != "") &&
            (r.I2 != "")
            ) return 2;
        else return 1;
    }
// Continuing validation, be sure the files exist and are not empty.
// If any file in a row is empty, make a note of it and proceed with the remainder

input_w_index_invalid_ch = Channel.create()
input_w_index_valid_ch = Channel.create()
input_w_index_ch.choice(
    input_w_index_invalid_ch,
    input_w_index_valid_ch
) {
    r -> (file(r.R1).isEmpty() || file(r.R2).isEmpty() || file(r.I1).isEmpty() || file(r.I2).isEmpty()) ? 0 : 1
}

input_no_index_invalid_ch = Channel.create()
input_no_index_valid_ch = Channel.create()
input_no_index_ch.choice(
    input_no_index_invalid_ch,
    input_no_index_valid_ch
) {
    r -> (file(r.R1).isEmpty() || file(r.R2).isEmpty()) ? 0 : 1
}

// For those with an index, make a channel for barcodecop
input_w_index_valid_ch
    .map{ sample -> [
        sample.specimen,
        file(sample.R1),
        file(sample.R2),
        file(sample.I1),
        file(sample.I2),
    ]}
    .set{ to_bcc_ch }
// ... and run barcodecop to validate the demultiplex
process barcodecop {
  container "golob/barcodecop:0.4.1__bcw_0.3.0"
  label 'io_limited'
  errorStrategy "retry"

  input:
    set specimen, file(R1), file(R2), file(I1), file(I2) from to_bcc_ch
  
  output:
    set specimen, file("${R1.getSimpleName()}_R1.bcc.fq.gz"), file("${R2.getSimpleName()}_R2.bcc.fq.gz") into post_bcc_ch
  """
  set -e

  barcodecop \
  ${I1} ${I2} \
  --match-filter \
  -f ${R1} \
  -o ${R1.getSimpleName()}_R1.bcc.fq.gz &&
  barcodecop \
  ${I1} ${I2} \
  --match-filter \
  -f ${R2} \
  -o ${R2.getSimpleName()}_R2.bcc.fq.gz
  """
}

// See if any read sets fail at this step
bcc_to_cutadapt_ch = Channel.create()
bcc_empty_ch = Channel.create()
post_bcc_ch
    .choice(
        bcc_to_cutadapt_ch,
        bcc_empty_ch
    ) {
        r -> (file(r[1]).isEmpty() || file(r[2]).isEmpty()) ? 1 : 0
    }

// Mix together the reads with no index with the reads with verified demultiplex

input_no_index_valid_ch
    .map { sample -> [
        sample.specimen,
        // Actual files
        file(sample.R1),
        file(sample.R2),
    ]}
    .mix(bcc_to_cutadapt_ch)
    .set{ demupltiplexed_ch}

// Step 2
process cutadapt {
  container "golob/cutadapt:2.3__bcw.0.3.0_al38B_FH"
  label 'io_limited'
  errorStrategy "retry"

  //publishDir "${params.output_folder}/noadapt/"

  input:
  set sample_name, file(R1), file(R2) from demupltiplexed_ch
  
  output:
  set sample_name, file("${R1.getSimpleName()}_R1.noadapt.fq.gz"), file("${R2.getSimpleName()}_R2.noadapt.fq.gz"), file("${R1.getSimpleName()}.cutadapt.log") into noadapt_ch

  """
  set -e 

  cutadapt \
  -j ${task.cpus} \
   -a ${params.adapter_F} -A ${params.adapter_R} \
  -o ${R1.getSimpleName()}_R1.noadapt.fq.gz -p ${R2.getSimpleName()}_R2.noadapt.fq.gz \
  ${R1} ${R1} > ${R1.getSimpleName()}.cutadapt.log
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
process removeHuman {
  container "golob/bwa:0.7.17__bcw.0.3.0I"
  errorStrategy "retry"
  label 'mem_veryhigh'


  input:
    file hg_index_tgz from hg_index_tgz
    set sample_name, file(R1), file(R2), file(cutadapt_log) from noadapt_ch
  
  output:
    set sample_name, file("${R1.getSimpleName()}.noadapt.nohuman.fq.gz"), file("${R2.getSimpleName()}.noadapt.nohuman.fq.gz"), file("${R1.getSimpleName()}.nohuman.log") into nohuman_ch

  afterScript "rm -rf hg_index/*"

  """
  set - e

  bwa_index_prefix=\$(tar -ztvf ${hg_index_tgz} | head -1 | sed \'s/.* //\' | sed \'s/.amb//\') && \
  echo BWA index file prefix is \${bwa_index_prefix} | tee -a ${R1.getSimpleName()}.nohuman.log && \
  echo Extracting BWA index | tee -a ${R1.getSimpleName()}.nohuman.log && \
  mkdir -p hg_index/ && \
  tar -I pigz -xf ${hg_index_tgz} -C hg_index/ | tee -a ${R1.getSimpleName()}.nohuman.log && \
  echo Files in index directory: | tee -a ${R1.getSimpleName()}.nohuman.log && \
  ls -l -h hg_index | tee -a ${R1.getSimpleName()}.nohuman.log && \
  echo Running BWA | tee -a ${R1.getSimpleName()}.nohuman.log && \
  bwa mem -t ${task.cpus} \
  -T ${params.min_hg_align_score} \
  -o alignment.sam \
  hg_index/\$bwa_index_prefix \
  ${R1} ${R2} \
  | tee -a ${R1.getSimpleName()}.nohuman.log && \
  echo Checking if alignment is empty  | tee -a ${R1.getSimpleName()}.nohuman.log && \
  [[ -s alignment.sam ]] && \
  echo Extracting Unaligned Pairs | tee -a ${R1.getSimpleName()}.nohuman.log && \
  samtools fastq alignment.sam \
  --threads ${task.cpus} -f 12 \
  -1 ${R1.getSimpleName()}.noadapt.nohuman.fq.gz -2 ${R2.getSimpleName()}.noadapt.nohuman.fq.gz \
  | tee -a ${R1.getSimpleName()}.nohuman.log && \
  echo Done | tee -a ${R1.getSimpleName()}.nohuman.log
  """
}

// Combine reads, grouping by specimen
nohuman_ch.groupTuple()
  .set{ sample_no_human_group_ch }

process combineReads {
    container = 'golob/fastatools:0.6.2__bcw.0.3.0'
    label = 'io_limited'

    publishDir path: "${params.output_folder}/qc/", mode: 'copy'

    input:
      set val(sample), file(R1s), file(R2s), file(logs) from sample_no_human_group_ch
    
    output:
      set val(sample), file("${sample}.R1.fastq.gz"), file("${sample}.R2.fastq.gz") into sample_qc_ch

    """
    combine_fastq_pairs.py \
    -1 ${R1s} \
    -2 ${R2s} \
    --normalize-ids \
    -o1 "${sample}.R1.fastq.gz" \
    -o2 "${sample}.R2.fastq.gz"
    """

}

// Output the final reads and manifest
sample_qc_ch.reduce('specimen,R1,R2\n'){ csvStr, row ->
            return  csvStr += "${row[0]},${params.output_folder}/qc/${row[1].name},${params.output_folder}/qc/${row[2].name}\n";
        }.set{manifestStr}

process outputManifest {
    container "golob/cutadapt:2.3__bcw.0.3.0_al38B_FH"

    publishDir path: "${params.output_folder}/", mode: 'copy'

    input:
        val manifestStr from manifestStr
    
    output:
        file 'manifest.qc.csv' into opManifest

    """
        echo "${manifestStr}" > manifest.qc.csv
    """
}
// */