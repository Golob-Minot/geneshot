#!/usr/bin/env nextflow
nextflow.enable.dsl=2

container__metaphlan2 = 'golob/metaphlan:2.8'
container_anndata = 'golob/python-anndata:0.9.2'

params.help = false
params.output = './results'
params.manifest = false

workflow Metaphlan2_wf {
    take:
        fastq_pairs
        fastq_unpaired
    
    main:
        // First handle the fastq / paired reads


        Metaphlan2_paired(fastq_pairs)
        Metaphlan2_unpaired(fastq_unpaired)
        Join_metaphlan2(
            Metaphlan2_paired.out.map{it[1]}.mix(
                Metaphlan2_unpaired.out.map{it[1]}
            ).toList()
        )
        MetaphlanLongToAnndata(
            Join_metaphlan2.out
        )
    
}

    
process Metaphlan2_paired {
    tag "MetaPhlAn2 composition for paired end fastq reads"
    container "${container__metaphlan2}"
    label = 'multithread'
    errorStrategy 'ignore'
    publishDir path: "${params.output_folder}MetaPhlAn2/by_specimen/", mode: "copy"

    input:
    tuple val(specimen), path(R1), path(R2)
    
    output:
    tuple val(specimen), path("${specimen}.metaphlan2.tsv")

"""
set -e

bowtie2 --threads ${task.cpus} -x /metaphlan/metaphlan_databases/mpa_v20_m200 -1 ${R1} -2 ${R2} | \
metaphlan2.py --nproc ${task.cpus} --input_type sam -t rel_ab_w_read_stats -o ${specimen}.metaphlan2.tsv

"""
}


process Metaphlan2_unpaired {
    tag "MetaPhlAn2 composition for unpaired fastq reads"
    container "${container__metaphlan2}"
    label = 'multithread'
    errorStrategy 'ignore'

    publishDir path: "${params.output_folder}MetaPhlAn2/by_specimen/", mode: "copy"

    input:
    tuple val(specimen), path(R1)
    
    output:
    tuple val(specimen), path("${specimen}.metaphlan2.tsv")

"""
set -e

bowtie2 --threads ${task.cpus} -x /metaphlan/metaphlan_databases/mpa_v20_m200 -f -U ${R1} | \
metaphlan2.py --nproc ${task.cpus} --input_type sam -t rel_ab_w_read_stats -o ${specimen}.metaphlan2.tsv

"""
}


process Join_metaphlan2 {
    container "${container_anndata}"
    label = 'io_mem'
    errorStrategy 'finish'
    publishDir path: "${params.output_folder}MetaPhlAn2/", mode: "copy"

    input:
    path metaphlan_tsv_list
    
    output:
    path "metaphlan2.results.csv.gz"

"""
#!/usr/bin/env python3

import pandas as pd

# Save all of the results to a dict
df = dict()

metaphlan_result_list = []
# Iterate over each input file
for fp in "${metaphlan_tsv_list.join(';;')}".split(";;"):

    # Make sure that the file has the expected suffix
    assert fp.endswith(".metaphlan2.tsv"), fp

    # Get the specimen name from the file name
    specimen_name = fp.replace(".metaphlan2.tsv", "")
    sp_df = pd.read_csv(
        fp,
        comment='#',
        names=[
            'clade_name',
            'relative_abundance',	
            'coverage',
            'average_genome_length_in_the_clade',
            'estimated_number_of_reads_from_the_clade'
        ],
        delimiter="\\t"
    )
    sp_df['specimen'] = specimen_name
    sp_df['clade_terminus'] = sp_df.clade_name.apply(lambda c: c.split("|")[-1])
    sp_df['rank'] = sp_df.clade_terminus.apply(lambda c: c.split('__')[0])
    try:
        sp_df['tax_name'] = sp_df.clade_terminus.apply(lambda c: c.split('__')[1])
    except IndexError:
        sp_df['tax_name'] = sp_df.clade_terminus.apply(lambda c: c.split('__')[0])

    metaphlan_result_list.append(sp_df)

pd.concat(metaphlan_result_list).to_csv("metaphlan2.results.csv.gz", index=None)
"""
    }

process MetaphlanLongToAnndata {
    container "${container_anndata}"
    label = 'io_mem'
    errorStrategy 'finish'
    publishDir path: "${params.output_folder}MetaPhlAn2/", mode: "copy"

    input:
    path "metaphlan2.results.csv.gz"
    
    output:
    path "metaphlan2.species.h5ad"

"""
#!/usr/bin/env python3

import pandas as pd
import anndata as ad
import numpy as np
import scipy

mpl = pd.read_csv('metaphlan2.results.csv.gz')
mpl_species = mpl[mpl['rank'] == 's']
all_specimens = sorted(set(mpl.specimen))
specimen_i = {sp: i for (i, sp) in enumerate(all_specimens)}
all_species = sorted(set(mpl[mpl['rank'] == 's'].tax_name))
species_i = {sp: i for (i, sp) in enumerate(all_species)}
mp_ad = ad.AnnData(
    scipy.sparse.coo_matrix(
        (
            mpl_species.estimated_number_of_reads_from_the_clade.values,
            (
                mpl_species.specimen.apply(specimen_i.get),
                mpl_species.tax_name.apply(species_i.get),
            )
        ),
        shape=(len(all_specimens), len(all_species)),
        dtype=np.int32
    ).tocsr(),
    obs=pd.DataFrame(index=all_specimens),
    var=pd.DataFrame(index=all_species),
)
mp_ad.write_h5ad("metaphlan2.species.h5ad")

"""
    }    

// Function which prints help message text
def helpMessage() {
    log.info"""
    A workflow oriented around obtaining the composition of a community via MetaPhlAn2.

    Usage:

    nextflow run modules/composition.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples (see below)

    Options:
      --output              Folder to place analysis outputs (default ./results)

      -w                    Working directory. Defaults to `./work`
    
    Manifest file:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. Each specimen should have exactly one row.
      Path to reads are specified by columns, `R1` and `R2`. These should be in fastq.gz format. 
      This workflow can handle single end reads. In that case only specify `R1` and omit or leave blank `R2` for single end.

    
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help || !params.manifest){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}



include { Read_manifest } from './general'
// Make sure that --output ends with trailing "/" characters
if (!params.output.endsWith("/")){
    params.output_folder = params.output.concat("/")
} else {
    params.output_folder = params.output
}


workflow {
    main:

    // Phase 0: Validation of input data
    manifest_file = Channel.from(file(params.manifest))
    // Read manifest splits out our manifest.
    manifest_qced = Read_manifest(manifest_file)
    

    // ################
    // # Composition  #
    // ################
    Metaphlan2_wf(
        manifest_qced.valid_paired.mix(
            manifest_qced.valid_paired_indexed
        ).map{
            r -> [r.specimen, file(r.R1), file(r.R2)]
        },
        manifest_qced.valid_unpaired.map{
            r -> [r.specimen, file(r.R1)]
        }
    )
    // */
}

