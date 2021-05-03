nextflow.enable.dsl=2

// Container versions
container__barcodecop = "quay.io/fhcrc-microbiome/barcodecop:barcodecop_0.5.3"
container__trimgalore = 'quay.io/biocontainers/trim-galore:0.6.6--0'
container__bwa = "quay.io/fhcrc-microbiome/bwa:bwa.0.7.17__bcw.0.3.0I"
container__fastatools = "quay.io/fhcrc-microbiome/fastatools:0.7.1__bcw.0.3.2"
container__ubuntu = "ubuntu:18.04"

// Function to filter a manifest to those rows which 
// have values for specimen, R1, and R2, but are missing any values in I1 or I2
def filter_no_index(manifest_ch){
    manifest_ch.filter { r ->
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
    }
}

// Function to filter a manifest to those rows which
// have values for specimen, R1, R2, I1, and I2
def filter_valid_index(manifest_ch){
    manifest_ch.filter { r ->
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
    }
}

workflow Preprocess_wf {
    take: indexed_ch
    take: paired_ch

    main:

    // Indexed rows into a channel for barcodecop
    indexed_ch.map{ sample -> [
        sample.specimen,
        file(sample.R1),
        file(sample.R2),
        file(sample.I1),
        file(sample.I2),
    ]}.set{ to_bcc_ch }

    // Run barcodecop
    Barcodecop(to_bcc_ch)

    // Raise an error if any of the samples fail barcodecop
    Barcodecop.out.bcc_empty_ch.filter {
        r -> (file(r[1]).isEmpty() || file(r[2]).isEmpty())
    }.map { r -> assert false: "Specimen failed barcodecop ${r.specimen}"}

    // Send the files which pass barcodecop to trim
    bcc_to_trim_ch = Barcodecop.out.bcc_to_cutadapt_ch
        .filter { 
                r -> (!file(r[1]).isEmpty() && !file(r[2]).isEmpty())
        }

    // Mix together the reads with no index with the reads with verified demultiplex
     paired_ch
        .map { sample -> [
            sample.specimen,
            // Actual files
            file(sample.R1),
            file(sample.R2),
        ]}
        .mix(bcc_to_trim_ch)
        .set{ demupltiplexed_ch}

    // Run catadapt
    TrimGalore(demupltiplexed_ch)

    // Download the human index if needed
    if (!params.hg_index) {
        Download_hg_index()
        hg_index_tgz = Download_hg_index.out
    } else {
        hg_index_tgz = file(params.hg_index)
    }

    // Remove the human reads
    BWA_remove_human(hg_index_tgz, TrimGalore.out)

    // Combine the reads by specimen name
    CombineReads(BWA_remove_human.out.groupTuple())

    // If the user specified --savereads, write out the manifest
    if (params.savereads) {
        WriteManifest(
            CombineReads.out
        )
    }

    // Set the outputs of the workflow
    emit:
        CombineReads.out

}




// Run barcodecop to validate the demultiplex
process Barcodecop {
    tag "Validate barcode demultiplexing for WGS reads"
    container "${container__barcodecop}"
    label 'multithread'
    errorStrategy 'finish'

    input:
        tuple val(specimen), file(R1), file(R2), file(I1), file(I2)

    output:
        tuple val(specimen), file("${R1}.bcc.fq.gz"), file("${R2}.bcc.fq.gz"), emit: bcc_to_cutadapt_ch
        tuple val(specimen), file("${R1}.bcc.fq.gz"), file("${R2}.bcc.fq.gz"), emit: bcc_empty_ch
"""
set -e

echo "Running barcodecop on ${R1}"
barcodecop \
${I1} ${I2} \
--match-filter \
-f ${R1} \
-o OUTPUT_R1.fq.gz
echo "Done"

echo "Renaming output file to ${R1}.bcc.fq.gz"
mv OUTPUT_R1.fq.gz ${R1}.bcc.fq.gz
echo "Done"

echo "Running barcodecop on ${R2}"
barcodecop \
${I1} ${I2} \
--match-filter \
-f ${R2} \
-o OUTPUT_R2.fq.gz
echo "Done"

echo "Renaming output file to ${R2}.bcc.fq.gz"
mv OUTPUT_R2.fq.gz ${R2}.bcc.fq.gz
echo "Done"
"""
}

// Use trim_galore to handle adapters / etc
process TrimGalore {
    container "${container__trimgalore}"
    label 'io_limited'
    errorStrategy 'ignore'

    input:
    tuple val(specimen), file(R1), file(R2)

    output:
    tuple val(specimen), file("${R1.getSimpleName()}__R1.tg.fastq.gz"), file("${R2.getSimpleName()}__R2.tg.fastq.gz")

    """
    set -e

    ln -s ${R1} R1.fastq.gz
    ln -s ${R2} R2.fastq.gz

    trim_galore \
    --gzip \
    --cores ${task.cpus} \
    --paired \
    R1.fastq.gz R2.fastq.gz

    mv R1_val_1.fq.gz "${R1.getSimpleName()}__R1.tg.fastq.gz"
    mv R2_val_2.fq.gz "${R2.getSimpleName()}__R2.tg.fastq.gz"
    """
}

// Process to download the human genome BWA index, already tarballed
process Download_hg_index {
    tag "Download human reference genome"
    container "${container__bwa}"
    errorStrategy "finish"
    label 'io_net'

    output:
    file 'hg_index.tar.gz'
    
"""
set -e

wget --quiet ${params.hg_index_url} -O hg_index.tar.gz
"""
}


// Process to remove human reads
process BWA_remove_human {
    tag "Remove human reads"
    container "${container__bwa}"
    errorStrategy 'finish'
    label 'mem_veryhigh'


    input:
        file hg_index_tgz
        tuple val(sample_name), file(R1), file(R2)

    output:
        tuple val(sample_name), file("${R1.getSimpleName()}.noadapt.nohuman.fq.gz"), file("${R2.getSimpleName()}.noadapt.nohuman.fq.gz")


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

rm -rf hg_index/*
echo Cleanup Done
"""
}

process JoinFASTQ {
    tag "Join FASTQ files per-specimen"
    container "${container__fastatools}"
    label = 'mem_medium'
    errorStrategy 'finish'

    // If the user sets --preprocess_output, write out the combined reads to that folder
    publishDir path: "${params.output}qc/", enabled: params.savereads, mode: "copy"

    input:
    tuple val(sample), file("R1.*.fastq.gz"), file("R2.*.fastq.gz")
    
    output:
    tuple val(sample), file("${sample}.R1.fastq.gz"), file("${sample}.R2.fastq.gz")

"""
set -e

ls -lah *

combine_fastq_pairs.py \
-1 R1*fastq.gz \
-2 R2*fastq.gz \
--normalize-ids \
-o1 "${sample}.R1.fastq.gz" \
-o2 "${sample}.R2.fastq.gz"

(( \$(gunzip -c "${sample}.R1.fastq.gz" | head | wc -l) > 1 ))
(( \$(gunzip -c "${sample}.R2.fastq.gz" | head | wc -l) > 1 ))

"""
}


process OutputManifest {
    container "${container__ubuntu}"

    publishDir path: "${params.output}qc/", enabled: params.savereads, mode: "copy"

    input:
        path R1s
        path R2s
        val manifestStr
    
    output:
        path 'manifest.qc.csv'
        path R1s
        path R2s

    """
        echo "${manifestStr}" > manifest.qc.csv
    """
}

// Workflow to publish a set of reads to a folder, along with a manifest
workflow WriteManifest {
    take:
        reads_ch

    main:
        // Make a manifest for the files in reads_ch
        // Output the final reads and manifest
        
        manifestStr = reads_ch.reduce(
            'specimen,R1,R2\n'
        ){ csvStr, row ->
            return  csvStr += "${row[0]},${params.output}qc/${row[1].name},${params.output}qc/${row[2].name}\n";
        }


        // Write the manifest CSV to a file
        OutputManifest(
            reads_ch.collect{ file(it[1]) },
            reads_ch.collect{ file(it[2]) },
            manifestStr
        )
        
}

workflow CombineReads {
    take:

        fastq_ch

    main:

        fastq_ch.branch {  // Split up the samples which have multiple FASTQ files
            single: it[1].size() == 1
            multiple: it[1].size() > 1
        }.set {
            grouped_fastq
        }

        JoinFASTQ(
            grouped_fastq.multiple
        )

    emit:
        grouped_fastq.single.map {
            r -> [r[0], r[1][0], r[2][0]]
        }.mix(
            JoinFASTQ.out
        )

}
//
// Steps to run preprocessing independently.
//

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below
params.nopreprocess = false
params.savereads = false
params.help = false
params.output = './results/'
params.manifest = null

// Preprocessing options
params.hg_index_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.bwa_index.tar.gz'
params.hg_index = false
params.min_hg_align_score = 30

// imports
include { Read_manifest } from './general'


// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run Golob-Minot/geneshot/preprocess <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    Options:
      --output              Folder to place analysis outputs (default ./results/)
      --savereads           If specified, save the preprocessed reads to the output folder (inside qc/)
      -w                    Working directory. Defaults to `./work`

    For preprocessing:
      --hg_index_url        URL for human genome index, defaults to current HG
      --hg_index            Cached copy of the bwa indexed human genome, TGZ format
      --min_hg_align_score  Minimum alignment score for human genome (default 30)

    Manifest:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This can be repeated. 
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      If index reads are provided, the column titles should be 'I1' and 'I2'

    """.stripIndent()
}


workflow {
    main:

 
    // Show help message if the user specifies the --help flag at runtime
    if (params.help || params.manifest == null){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 0
    }
    // Read and validate manifest
    manifest_file = Channel.from(file(params.manifest))
    manifest_qced = Read_manifest(manifest_file)
    // Actually preprocess

    Preprocess_wf(
        manifest_qced.valid_paired_indexed,
        manifest_qced.valid_paired
    )

}