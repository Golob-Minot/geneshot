// Container versions
container__barcodecop = "quay.io/fhcrc-microbiome/barcodecop:barcodecop_0.5.3"
container__cutadapt = "quay.io/fhcrc-microbiome/cutadapt:cutadapt_2.3_bcw_0.3.1"
container__bwa = "quay.io/fhcrc-microbiome/bwa:bwa.0.7.17__bcw.0.3.0I"

// Input to this workflow is a manifest CSV 



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
     paired_ch
        .map { sample -> [
            sample.specimen,
            // Actual files
            file(sample.R1),
            file(sample.R2),
        ]}
        .mix(bcc_to_cutadapt_ch)
        .set{ demupltiplexed_ch}

    // Run catadapt
    cutadapt(demupltiplexed_ch)

    // Download the human index if needed
    if (!params.hg_index) {
        download_hg_index()
        hg_index_tgz = download_hg_index.out
    } else {
        hg_index_tgz = file(params.hg_index)
    }

    // Remove the human reads
    bwa(hg_index_tgz, cutadapt.out)

    // Set the outputs of the workflow
    emit:
        bwa.out

}




// Run barcodecop to validate the demultiplex
process barcodecop {
    tag "Validate barcode demultiplexing for WGS reads"
    container "${container__barcodecop}"
    label 'mem_medium'
    maxRetries 10

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

// Process to run catadapt
process cutadapt {
    tag "Trim adapters from WGS reads"
    container "${container__cutadapt}"
    label 'mem_medium'
    maxRetries 10

    input:
    tuple val(sample_name), file(R1), file(R2)

    output:
    tuple val(sample_name), file("${R1.getSimpleName()}_R1.noadapt.fq.gz"), file("${R2.getSimpleName()}_R2.noadapt.fq.gz"), file("${R1.getSimpleName()}.cutadapt.log")

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
    tag "Download human reference genome"
    container "${container__bwa}"
    label 'io_limited'

    output:
    file 'hg_index.tar.gz'
    
"""#!/bin/bash
set -e

wget --quiet ${params.hg_index_url} -O hg_index.tar.gz

# Make sure that the entire file was downloaded
tar -tzvf hg_index.tar.gz
"""
}


// Process to remove human reads
process bwa {
    tag "Remove human reads"
    container "${container__bwa}"
    maxRetries 10
    label 'mem_veryhigh'


    input:
        file hg_index_tgz
        tuple val(sample_name), file(R1), file(R2), file(cutadapt_log)

    output:
        tuple val(sample_name), file("${R1.getSimpleName()}.noadapt.nohuman.fq.gz"), file("${R2.getSimpleName()}.noadapt.nohuman.fq.gz")

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