container__metaphlan2 = 'golob/metaphlan:2.8'


workflow composition_wf {
    take:
        fastq_pairs
        fasta_files
    
    main:
        // First handle the fastq / paired reads


        metaphlan2_fastq(fastq_pairs)
        metaphlan2_fasta(fasta_files)

    
}

    
    process metaphlan2_fastq {
    tag "MetaPhlAn2 composition for paired end fastq reads"
    container "${container__metaphlan2}"
    label = 'mem_veryhigh'
    //errorStrategy 'retry'
    //maxRetries 10

    // If the user sets --preprocess_output, write out the combined reads to that folder
    publishDir path: "${params.output_folder}MetaPhlAn2/", mode: "copy"

    input:
    tuple specimen, file(R1), file(R2)
    
    output:
    tuple val(specimen), file("${specimen}.metaphlan2.tsv")

"""
set -e

bowtie2 --threads ${task.cpus} -x /metaphlan/metaphlan_databases/mpa_v20_m200 -1 ${R1} -2 ${R2} | \
metaphlan2.py --nproc ${task.cpus} --input_type sam -t rel_ab_w_read_stats -o ${specimen}.metaphlan2.tsv

"""
    }


    process metaphlan2_fasta {
    tag "MetaPhlAn2 composition for paired end fastq reads"
    container "${container__metaphlan2}"
    label = 'mem_veryhigh'
    //errorStrategy 'retry'
    //maxRetries 10

    // If the user sets --preprocess_output, write out the combined reads to that folder
    publishDir path: "${params.output_folder}MetaPhlAn2/", mode: "copy"

    input:
    tuple specimen, file(R1)
    
    output:
    tuple val(specimen), file("${specimen}.metaphlan2.tsv")

"""
set -e

bowtie2 --threads ${task.cpus} -x /metaphlan/metaphlan_databases/mpa_v20_m200 -f -U ${R1} | \
metaphlan2.py --nproc ${task.cpus} --input_type sam -t rel_ab_w_read_stats -o ${specimen}.metaphlan2.tsv

"""
    }