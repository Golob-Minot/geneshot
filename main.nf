#!/usr/bin/env nextflow

/*
  Geneshot: A pipeline to robustly identify which alleles (I.e. peptide coding sequences)
  are present in a microbial community.

  There is an optional preprocessing submodule, 
  followed by assembly (spades) + extraction of peptide coding sequences from the contigs
  The short reads are then aligned against the assembled peptides plus uniref100.
  We use the FAMLI algorithm to adjuticate these alignments.
  Annotations can follow. 

  I. (Optional) Geneshot preprocessing submodule:
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
params.nopreprocess = false
params.help = false
params.output = './results'

// Preprocessing options
params.adapter_F = "CTGTCTCTTATACACATCT"
params.adapter_R = "CTGTCTCTTATACACATCT"
params.hg_index_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30

// Assembly options
params.phred_offset = 33 // spades
params.centre = 'geneshot' // prokka

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run jgolob/geneshot <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    Options:
      --output              Folder to place outputs (default ./results)
      -w                    Working directory. Defaults to `./work`
    For preprocessing:
      --hg_index_url        URL for human genome index, defaults to current HG
      --hg_index            Cached copy of the bwa indexed human genome, TGZ format
      --adapter_F           Forward sequencing adapter sequence (to be removed)
      --adapter_R           Reverse sequencing adapter sequence (to be removed)
                              (Adapter sequences default to nextera adapters)
      --min_hg_align_score  Minimum alignment score for human genome (default 30)
    For Assembly:
      --phred_offset        for spades. Default 33.
      --centre              Centre for use in prokka. default = 'geneshot'
      

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

// Phase I: Preprocessing
if (!params.nopreprocess) {
    // Implicit else we have a manifest. Time to do some validation:
    // Split into index provided or not. Further validate critical elements are present
    // e.g. A specimen, R1 and R2

    Channel.from(file(params.manifest))
        .splitCsv(header: true, sep: ",")
        .into {
            input_w_index_ch;
            input_no_index_ch;
            input_invalid_ch
        }

    input_invalid_ch
        .filter { r ->
            (r.specimen == null) ||
            (r.R1 == null) ||
            (r.R2 == null) ||
            (r.specimen == "") ||
            (r.R1 == "") ||
            (r.R2 == "")        
        }.set { input_invalid_ch }

    input_no_index_ch
        .filter { r ->
            (r.specimen != null) &&
            (r.R1 != null) &&
            (r.R2 != null) &&
            (r.specimen != "") &&
            (r.R1 != "") &&
            (r.R2 != "") &&
            (
                (r.I1 == null) ||
                (r.I2 == null) ||
                (r.I1 == "") ||
                (r.I2 == "")
            )
        }.into { 
            input_no_index_valid_ch;
            input_no_index_invalid_ch;
        }

    // be sure files exist, are not empty.
    input_no_index_valid_ch
        .filter {
            r -> (!file(r.R1).isEmpty() && !file(r.R2).isEmpty())
        }.set {
            input_no_index_valid_ch
        }
    input_no_index_invalid_ch
        .filter {
            r -> (file(r.R1).isEmpty() || file(r.R2).isEmpty())
        }.set {
            input_no_index_invalid_ch
        }

    input_w_index_ch
        .filter { r ->
            (r.specimen != null) &&
            (r.R1 != null) &&
            (r.R2 != null) &&
            (r.specimen != "") &&
            (r.R1 != "") &&
            (r.R2 != "") &&
            (r.I1 != null) &&
            (r.I2 != null) &&
            (r.I1 != "") &&
            (r.I2 != "")
        }.into { 
            input_w_index_invalid_ch;
            input_w_index_valid_ch;
        }
    // be sure files exist, are not empty.
    input_w_index_valid_ch
        .filter {
            r -> (!file(r.R1).isEmpty() && !file(r.R2).isEmpty() && !file(r.I1).isEmpty() && !file(r.I2).isEmpty())
        }.set {
            input_w_index_valid_ch
        }
    input_w_index_invalid_ch
        .filter {
            r -> (file(r.R1).isEmpty() || file(r.R2).isEmpty() || file(r.I1).isEmpty() || file(r.I2).isEmpty())
        }.set {
            input_w_index_invalid_ch
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
    errorStrategy 'retry'
    maxRetries 10

    input:
        set specimen, file(R1), file(R2), file(I1), file(I2) from to_bcc_ch
    
    output:
        set specimen, file("${R1.getSimpleName()}_R1.bcc.fq.gz"), file("${R2.getSimpleName()}_R2.bcc.fq.gz") into bcc_to_cutadapt_ch
        set specimen, file("${R1.getSimpleName()}_R1.bcc.fq.gz"), file("${R2.getSimpleName()}_R2.bcc.fq.gz") into bcc_empty_ch
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
    bcc_to_cutadapt_ch
        .filter { 
             r -> (!file(r[1]).isEmpty() && !file(r[2]).isEmpty())
        }
        .set { bcc_to_cutadapt_ch }
    bcc_empty_ch
        .filter {
            r -> (file(r[1]).isEmpty() || file(r[2]).isEmpty())
        }
        .set { bcc_empty_ch }
    // Mix together the reads with no index with the reads with verified demultiplex
    input_no_index_valid_ch
        .map { sample -> [
            sample.specimen,
            // Actual files
            file(sample.R1),
            file(sample.R2),
        ]}
        .mix(bcc_to_cutadapt_ch)
        .set{ demupltiplexed_ch }

    // Step 2
    process cutadapt {
    container "golob/cutadapt:2.3__bcw.0.3.0_al38B_FH"
    label 'io_limited'
    errorStrategy 'retry'
    maxRetries 10

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
    errorStrategy 'retry'
    maxRetries 10
    label 'mem_veryhigh'


    input:
        file hg_index_tgz from hg_index_tgz
        set sample_name, file(R1), file(R2), file(cutadapt_log) from noadapt_ch
    
    output:
        set sample_name, file("${R1.getSimpleName()}.noadapt.nohuman.fq.gz"), file("${R2.getSimpleName()}.noadapt.nohuman.fq.gz"), file("${R1.getSimpleName()}.nohuman.log") into nohuman_ch

    afterScript "rm -rf hg_index/*"

    """
    set - e


    bwa_index_fn=\$(tar -ztvf ${hg_index_tgz} | head -1 | sed \'s/.* //\')
    bwa_index_prefix=\${bwa_index_fn%.*}
    echo BWA index prefix is \${bwa_index_prefix} | tee -a ${R1.getSimpleName()}.nohuman.log
    echo Extracting BWA index | tee -a ${R1.getSimpleName()}.nohuman.log
    mkdir -p hg_index/ 
    tar -I pigz -xf ${hg_index_tgz} -C hg_index/ | tee -a ${R1.getSimpleName()}.nohuman.log
    echo Files in index directory: | tee -a ${R1.getSimpleName()}.nohuman.log
    ls -l -h hg_index | tee -a ${R1.getSimpleName()}.nohuman.log
    echo Running BWA | tee -a ${R1.getSimpleName()}.nohuman.log
    bwa mem -t ${task.cpus} \
    -T ${params.min_hg_align_score} \
    -o alignment.sam \
    hg_index/\$bwa_index_prefix \
    ${R1} ${R2} \
    | tee -a ${R1.getSimpleName()}.nohuman.log
    echo Checking if alignment is empty  | tee -a ${R1.getSimpleName()}.nohuman.log
    [[ -s alignment.sam ]]
    echo Extracting Unaligned Pairs | tee -a ${R1.getSimpleName()}.nohuman.log
    samtools fastq alignment.sam \
    --threads ${task.cpus} -f 12 \
    -1 ${R1.getSimpleName()}.noadapt.nohuman.fq.gz -2 ${R2.getSimpleName()}.noadapt.nohuman.fq.gz \
    | tee -a ${R1.getSimpleName()}.nohuman.log
    echo Done | tee -a ${R1.getSimpleName()}.nohuman.log
    """
    }

    // Combine reads, grouping by specimen
    nohuman_ch.groupTuple()
    .set{ sample_no_human_group_ch }

    process combineReads {
        container = 'golob/fastatools:0.7.0__bcw.0.3.0'
        label = 'io_limited'
        errorStrategy 'retry'
        maxRetries 10

        publishDir path: "${params.output}/qc/", mode: 'copy'

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
    sample_qc_ch
        .into {
            to_qc_manifest_ch
            to_assembly_ch
            to_interleave_ch
        }

    // Output the final reads and manifest
    to_qc_manifest_ch.reduce('specimen,R1,R2\n'){ csvStr, row ->
                return  csvStr += "${row[0]},${params.output}/qc/${row[1].name},${params.output}/qc/${row[2].name}\n";
            }.set{manifestStr}

    process outputManifest {
        container "golob/cutadapt:2.3__bcw.0.3.0_al38B_FH"

        publishDir path: "${params.output}/", mode: 'copy'

        input:
            val manifestStr from manifestStr
        
        output:
            file 'manifest.qc.csv' into opManifest

        """
            echo "${manifestStr}" > manifest.qc.csv
        """
    }
} // End optional phase I: Preprocessing
// TODO an else to read in the manifest for assembly, remaining if one is not doing preprocessing

// Many of the processes from here request the paired reads to be interleaved. We will do just that
process interleave {
    container "golob/fastatools:0.7.0__bcw.0.3.0"
    label 'io_limited'
    errorStrategy "retry"

    input:
    set sample_name, file(fastq1), file(fastq2) from to_interleave_ch

    output:
    set sample_name, file("${fastq1.getSimpleName()}.interleaved.fastq.gz") into concatenated_ch

    """
    set -e

    # Some basic checks that the files exist and the line numbers match
    [[ -s "${fastq1}" ]]
    [[ -s "${fastq2}" ]]
    (( \$(gunzip -c ${fastq1} | wc -l) == \$(gunzip -c ${fastq2} | wc -l) ))

    # Now interleave the files
    paste <(gunzip -c ${fastq1}) <(gunzip -c ${fastq2}) | paste - - - - | awk -v OFS="\\n" -v FS="\\t" '{print(\$1,\$3,\$5,\$7,\$2,\$4,\$6,\$8)}' | gzip -c > "${fastq1.getSimpleName()}.interleaved.fastq.gz"
    """
}
concatenated_ch
    .into {
        to_count_ch;
        to_metaphlan_ch
    }

// Count the number of input reads
process countReads {
  container "golob/fastatools:0.7.0__bcw.0.3.0"
  cpus 1
  memory "4 GB"
  errorStrategy "retry"
  
  input:
  set sample_name, file(fastq) from to_count_ch
  
  output:
  file "${sample_name}.countReads.csv" into total_counts

  """
set -e

[[ -s ${fastq} ]]

n=\$(gunzip -c "${fastq}" | awk 'NR % 4 == 1' | wc -l)
echo "${sample_name},\$n" > "${sample_name}.countReads.csv"
  """
}
// Make a single file which summarizes the number of reads across all samples
// This is only run after all of the samples are done processing through the
// 'total_counts' channel, which is transformed by the .collect() command into
// a single list containing all of the data from all samples.
process countReadsSummary {
    container "golob/fastatools:0.7.0__bcw.0.3.0"
    // The output from this process will be copied to the --output_folder specified by the user
    publishDir "${params.output}/qc/", mode: 'copy'
    errorStrategy "retry"

    input:
    // Because the input channel has been collected into a single list, this process will only be run once
    file readcount_csv_list from total_counts.collect()

    output:
    file "readcounts.csv" into readcounts_csv


    """
    set -e

    echo name,n_reads > readcounts.csv
    cat ${readcount_csv_list} >> readcounts.csv
    """
}
// PHASE II: Basic details of the reads / communities
// qc read count by specimen
// metaphlan
process metaphlan2 {
    container "quay.io/fhcrc-microbiome/metaphlan@sha256:51b416458088e83d0bd8d840a5a74fb75066b2435d189c5e9036277d2409d7ea"
    label = 'multithread'
    errorStrategy 'finish'
    maxRetries 10
    publishDir "${params.output}/metaphlan2/${sample_name}/", mode: 'copy'

    input:
    set val(sample_name), file(input_fastq) from to_metaphlan_ch
    
    output:
    file "${sample_name}.metaphlan.tsv" into metaphlan_for_summary, metaphlan_for_humann

    """
    set -e
    metaphlan2.py \
    --input_type fastq --tmp_dir ./ -o ${sample_name}.metaphlan.tsv ${input_fastq}
    """
}

// PHASE III: Assembly!
// Assembly with metaspades
process metaspadesAssembly {
    container 'golob/spades:3.13.1__bcw.0.3.1'
    label 'mem_veryhigh'
    //errorStrategy "retry"

    publishDir "${params.output}/assembly/${specimen}", mode: 'copy'

    input:
        set specimen, file(R1), file(R2) from to_assembly_ch
    
    output:
        set specimen, file("${specimen}.contigs.fasta"), file("${specimen}.scaffolds.fasta"), file("${specimen}.metaspades.log") into assembly_ch
    
    """
    set -e 
    echo 
    metaspades.py \
    --meta \
    --phred-offset ${params.phred_offset} \
    -1 ${R1} -2 ${R2} \
    -o . \
    -t ${task.cpus} -m ${task.memory.toMega() / 1024} | 
    tee -a ${specimen}.metaspades.log
    mv contigs.fasta ${specimen}.contigs.fasta
    mv scaffolds.fasta ${specimen}.scaffolds.fasta
    """
}

// Annotation with prokka

process prokkaAnnotate {
    container 'golob/prokka:1.1.14__bcw.0.3.1'
    label 'mem_veryhigh'
    errorStrategy "retry"
    publishDir "${params.output}/prokka/${specimen}/", mode: 'copy'

    input:
        set val(specimen) , file(contigs), file(scaffolds), file(spades_log) from assembly_ch
    
    output:
        set \
            val(specimen), \
            file("prokka/${specimen}.tbl.gz"), \
            file("prokka/${specimen}.err.gz"), \
            file("prokka/${specimen}.faa.gz"), \
            file("prokka/${specimen}.ffn.gz"), \
            file("prokka/${specimen}.fsa.gz"), \
            file("prokka/${specimen}.gff.gz"), \
            file("prokka/${specimen}.log.gz"), \
            file("prokka/${specimen}.sqn.gz"), \
            file("prokka/${specimen}.tsv.gz"), \
            file("prokka/${specimen}.txt.gz") into prokka_ch

        set val(specimen), file("prokka/${specimen}.faa.gz") into specimen_prokka_faa_ch
    
    """
    set -e 

    prokka \
    --outdir prokka/ \
    --prefix ${specimen} \
    --cpus ${task.cpus} \
    --metagenome \
    --compliant \
    --centre ${params.centre} \
    ${scaffolds} &&
    ls -l -h prokka/ &&
    gzip prokka/${specimen}.tbl &&
    gzip prokka/${specimen}.err &&
    gzip prokka/${specimen}.faa &&
    gzip prokka/${specimen}.ffn &&
    gzip prokka/${specimen}.fna &&
    gzip prokka/${specimen}.fsa &&
    gzip prokka/${specimen}.gff &&
    gzip prokka/${specimen}.log &&
    gzip prokka/${specimen}.sqn &&
    gzip prokka/${specimen}.tsv &&
    gzip prokka/${specimen}.txt
    """
}

// PHASE IV: Alignment
/*
today = new Date().format("yyyy-MM-dd")

// First download the uniref100 library
process dlUniref100 {
    container "golob/fastatools:0.7.0__bcw.0.3.0"
    label "io_limited"
    input:
        val today

    output:
        file "uniref100__${today}.fasta.gz" into uniref100_fasta_f
    """
    set -e
    
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
    mv uniref100.fasta.gz uniref100__${today}.fasta.gz
    """
}

process diamondIndexUR100 {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label "mem_veryhigh"

}

// Diamond
// FAMLI



// */