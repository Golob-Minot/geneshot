
process combineReads {
    container = 'golob/fastatools:0.7.0__bcw.0.3.0'
    label = 'io_limited'
    // errorStrategy 'retry'
    maxRetries 10

    // If the user sets --preprocess_output, write out the combined reads to that folder
    publishDir path: "${params.preprocess_output_folder}", enabled: params.preprocess_output_folder != false

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
"""

}

process outputManifest {
    container "golob/cutadapt:2.3__bcw.0.3.0_al38B_FH"

    publishDir path: "${params.preprocess_output_folder}"

    input:
        val manifestStr
    
    output:
        file 'manifest.qc.csv'

    """
        echo "${manifestStr}" > manifest.qc.csv
    """
}

// Workflow to publish a set of reads to a folder, along with a manifest
workflow writeManifest {
    get:
        reads_ch

    main:
        // Make a manifest for the files in reads_ch
        // Output the final reads and manifest

        
        manifestStr = reads_ch.reduce(
            'specimen,R1,R2\n'
        ){ csvStr, row ->
            return  csvStr += "${row[0]},${params.preprocess_output_folder}${row[1].name},${params.preprocess_output_folder}${row[2].name}\n";
        }

        // Write the manifest CSV to a file
        outputManifest(manifestStr)
        
}