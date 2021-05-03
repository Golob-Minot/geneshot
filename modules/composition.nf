container__metaphlan2 = 'golob/metaphlan:2.8'
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"


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
    label = 'mem_medium'
    errorStrategy 'finish'

    // If the user sets --preprocess_output, write out the combined reads to that folder
    publishDir path: "${params.output_folder}MetaPhlAn2/", mode: "copy"

    input:
    tuple val(specimen), file(R1), file(R2)
    
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
    errorStrategy 'finish'
    //maxRetries 10

    // If the user sets --preprocess_output, write out the combined reads to that folder
    publishDir path: "${params.output_folder}MetaPhlAn2/", mode: "copy"

    input:
    tuple val(specimen), file(R1)
    
    output:
    tuple val(specimen), file("${specimen}.metaphlan2.tsv")

"""
set -e

bowtie2 --threads ${task.cpus} -x /metaphlan/metaphlan_databases/mpa_v20_m200 -f -U ${R1} | \
metaphlan2.py --nproc ${task.cpus} --input_type sam -t rel_ab_w_read_stats -o ${specimen}.metaphlan2.tsv

"""
    }


    process join_metaphlan2 {
    container "${container__pandas}"
    label = 'mem_medium'
    errorStrategy 'finish'

    input:
    path metaphlan_tsv_list
    
    output:
    path "metaphlan.results.csv.gz"

"""
#!/usr/bin/env python3

import pandas as pd

# Save all of the results to a dict
df = dict()

# Iterate over each input file
for fp in "${metaphlan_tsv_list}".split(" "):

    # Make sure that the file has the expected suffix
    assert fp.endswith(".metaphlan2.tsv"), fp

    # Get the specimen name from the file name
    specimen_name = fp.replace(".metaphlan2.tsv", "")

    df[specimen_name] = pd.read_csv(
        fp,
        skiprows=1,
        sep="\\t"
    ).set_index(
        "#clade_name"
    )[
        "relative_abundance"
    ]

# Format all results as a DataFrame
df = pd.DataFrame(df)

# Fill in all missing values with 0
df = df.fillna(0)

# Save to a single CSV
df.reset_index(
).rename(
    columns = dict([("index", "taxon")])
).to_csv("metaphlan.results.csv.gz")
    

"""
    }
