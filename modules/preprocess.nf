// Input to this workflow is a manifest CSV 

// Function to read in a CSV and return a Channel
def read_manifest(manifest_file){
    manifest_file.splitCsv(
        header: true, 
        sep: ","
    )
}

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

workflow preprocess_wf {

    get:
    manifest_file

    main:

    // Start by checking the manifest for rows which are invalid in any way
    
    // Raise an error if there are any rows which are missing any values for specimen, R1, or R2
    read_manifest(manifest_file).filter { r ->
        (r.specimen == null) ||
        (r.R1 == null) ||
        (r.R2 == null) ||
        (r.specimen == "") ||
        (r.R1 == "") ||
        (r.R2 == "")        
    }.count(
    ).map { assert it == 0: "Found lines missing values for specimen, R1, or R2"}

    // Raise an error if there are any rows which point to empty or missing R1 / R2 files
    filter_no_index(
        read_manifest(manifest_file)
    ).filter {
        r -> (file(r.R1).isEmpty() || file(r.R2).isEmpty())
    }.count(
    ).map { assert it == 0: "Found lines pointing to empty files for R1 or R2"}

    // Raise an error if there are any rows which point to empty or missing files for R1, R2, I1, or I2
    filter_valid_index(
        read_manifest(manifest_file)
    ).filter {
        r -> (file(r.R1).isEmpty() || file(r.R2).isEmpty() || file(r.I1).isEmpty() || file(r.I2).isEmpty())
    }.count(
    ).map { assert it == 0: "Found lines pointing to empty files for R1, R2, I1, or I2"}

    // Get rows which have valid R1 and R2, but no values provided for I1 or I2
    input_no_index_valid_ch = filter_no_index(
        read_manifest(manifest_file)
    ).filter {
        r -> (!file(r.R1).isEmpty() && !file(r.R2).isEmpty())
    }

    // Get rows which have valid R1, R2, I1, and I2
    // and make a channel for barcodecop
    to_bcc_ch = filter_valid_index(
        read_manifest(manifest_file)
    ).filter {
        r -> (!file(r.R1).isEmpty() && !file(r.R2).isEmpty() && !file(r.I1).isEmpty() && !file(r.I2).isEmpty())
    }.map{ sample -> [
        sample.specimen,
        file(sample.R1),
        file(sample.R2),
        file(sample.I1),
        file(sample.I2),
    ]}

    // Run barcodecop
    barcodecop(to_bcc_ch)

    // Raise an error if any of the samples fail barcodecop
    barcodecop.out.bcc_empty_ch.filter {
        r -> (file(r[1]).isEmpty() || file(r[2]).isEmpty())
    }.map { r -> assert false: "Specimen failed barcodecop ${r.specimen}"}

    // Send the files which pass barcodecop to cutadapt
    bcc_to_cutadapt_ch = barcodecop.out.bcc_to_cutadapt_ch
        .filter { 
                r -> (!file(r[1]).isEmpty() && !file(r[2]).isEmpty())
        }

    // Mix together the reads with no index with the reads with verified demultiplex
    demupltiplexed_ch = input_no_index_valid_ch
        .map { sample -> [
            sample.specimen,
            // Actual files
            file(sample.R1),
            file(sample.R2),
        ]}
        .mix(bcc_to_cutadapt_ch)

    // Run catadapt
    cutadapt(demupltiplexed_ch)

    // Download the human index if needed
    if (!params.hg_index) {
        hg_index_tgz = download_hg_index().out
    } else {
        hg_index_tgz = file(params.hg_index)
    }

    // Remove the human reads
    removeHuman(hg_index_tgz, cutadapt.out)

    // Set the outputs of the workflow
    emit:
        removeHuman.out

}




// Run barcodecop to validate the demultiplex
process barcodecop {
    container "golob/barcodecop:0.4.1__bcw_0.3.0"
    label 'io_limited'
    errorStrategy 'retry'
    maxRetries 10

    input:
        tuple specimen, file(R1), file(R2), file(I1), file(I2)

    output:
        tuple specimen, file("${R1.getSimpleName()}_R1.bcc.fq.gz"), file("${R2.getSimpleName()}_R2.bcc.fq.gz"), emit: bcc_to_cutadapt_ch
        tuple specimen, file("${R1.getSimpleName()}_R1.bcc.fq.gz"), file("${R2.getSimpleName()}_R2.bcc.fq.gz"), emit: bcc_empty_ch
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

// Process to run catadapt
process cutadapt {
    container "golob/cutadapt:2.3__bcw.0.3.0_al38B_FH"
    label 'mem_medium'
    errorStrategy 'retry'
    maxRetries 10

    input:
    tuple sample_name, file(R1), file(R2)

    output:
    tuple sample_name, file("${R1.getSimpleName()}_R1.noadapt.fq.gz"), file("${R2.getSimpleName()}_R2.noadapt.fq.gz"), file("${R1.getSimpleName()}.cutadapt.log")

"""
set -e 

cutadapt \
-j ${task.cpus} \
-a ${params.adapter_F} -A ${params.adapter_R} \
-o ${R1.getSimpleName()}_R1.noadapt.fq.gz -p ${R2.getSimpleName()}_R2.noadapt.fq.gz \
${R1} ${R1} > ${R1.getSimpleName()}.cutadapt.log
"""
}

// Process to download the human genome BWA index, already tarballed
process download_hg_index {
    container "golob/bwa:0.7.17__bcw.0.3.0I"
    errorStrategy "retry"
    label 'io_limited'

    output:
    file 'hg_index.tar.gz'
    
"""
set -e

wget ${params.hg_index_url} -O hg_index.tar.gz
"""
}


// Process to remove human reads
process removeHuman {
    container "golob/bwa:0.7.17__bcw.0.3.0I"
    errorStrategy 'retry'
    maxRetries 10
    label 'mem_veryhigh'


    input:
        file hg_index_tgz
        tuple sample_name, file(R1), file(R2), file(cutadapt_log)

    output:
        tuple sample_name, file("${R1.getSimpleName()}.noadapt.nohuman.fq.gz"), file("${R2.getSimpleName()}.noadapt.nohuman.fq.gz")

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