#!/usr/bin/env nextflow

/*
  Geneshot assembly submodule:
    Steps:
    1) Assemble into contigs with metaspades
    2) Prokka to extract features
    3) Eggnogmapper to identify features
*/

// Default values for boolean flags
// If these are not set by the user, then they will be set to the values below
// This is useful for the if/then control syntax below

params.output_folder = '.'
params.help = false

// spades
params.phred_offset = 33

//
params.centre = 'geneshot'

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run geneshot-assembly.nf <ARGUMENTS>
    
    Required Arguments:
      --manifest            CSV file listing samples to preprocess (see below)

    Options:
      --output_folder       Folder to place outputs (default invocation dir)
      -w                    Working directory. Defaults to `./work`
      --phred_offset        for spades. Default 33.
      --centre              Centre for use in prokka. default = 'geneshot'


    Batchfile:
      The manifest is a CSV with a header indicating which samples correspond to which files.
      The file must contain a column `specimen`. This can be repeated. 
      Data is only accepted as paired reads.
      Reads are specified by columns, `R1` and `R2`.
      The presumption here is these reads are QCed and reads for a specimen are pre-combined
      into one pair. Each specimen should have only one row.
    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
if (params.help){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 0
}

Channel.from(file(params.manifest))
        .splitCsv(header: true, sep: ",")
        .map { sample -> [
          sample.specimen,
          // Actual files
          file(sample.R1),
          file(sample.R2),
        ]}
        .set{ input_ch }

// Assembly with metaspades
process metaspadesAssembly {
    container 'golob/spades:3.13.1__bcw.0.3.1'
    label 'mem_veryhigh'
    //errorStrategy "retry"

    publishDir "${params.output_folder}/${specimen}/assembly/", mode: 'copy'

    input:
        set specimen, file(R1), file(R2) from input_ch
    
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
    publishDir "${params.output_folder}/${specimen}/prokka/", mode: 'copy'

    input:
        set val(specimen), file(contigs), file(scaffolds), file(spades_log) from assembly_ch
    
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

/*
process eggnogMapperDownloadDB__Diamond {
    container 'golob/eggnog-mapper:1.0.3__bcw.0.3.1B'
    //label = 'multithread'
    label = 'io_limited'
    errorStrategy "ignore"
    afterScript "rm -rf /tmp/enmdb/"

    /*
    output:
        file("eggnog.db.gz") into eggnog_db_f
        file("eggnog_proteins.dmnd.gz") into eggnog_proteins_dmnd_f  
    """
    set -e
    
    mkdir -p /tmp/enmdb/ &&
    download_eggnog_data.py none -y --data_dir /tmp/enmdb/ &&
    pigz -p ${task.cpus} -c /tmp/db/emdb/eggnog.db > eggnog.db.gz
    pigz -p ${task.cpus} -c /tmp/db/emdb/eggnog_proteins.dmnd > eggnog_proteins.dmnd.gz
    """

    output:
        file("eggnog.db") into eggnog_db_f
        file("eggnog_proteins.dmnd") into eggnog_proteins_dmnd_f  

    """
    set -e
    
    mkdir -p /tmp/enmdb/ &&
    download_eggnog_data.py none -y --data_dir /tmp/enmdb/ &&
    ls -l -h /tmp/db/emdb && 
    tree /tmp/db/ &&
    mv /tmp/db/emdb/eggnog.db eggnog.db &&
    mv /tmp/db/emdb/eggnog_proteins.dmnd eggnog_proteins.dmnd
    """
}

process eggnogMap {
    container 'golob/eggnog-mapper:1.0.3__bcw.0.3.1B'
    label 'mem_veryhigh'
    errorStrategy "retry"

    input:
        set val(specimen), file(faa) from specimen_prokka_faa_ch
        file(eggnog_db_f)
        file(eggnog_proteins_dmnd_f)

    output:
        set val(specimen), file(" ${specimen}.egm.emapper.annotations.gz") into eggnogmap_ch

    """
    set -e 

    emapper.py \
    -i ${faa} \
    -m diamond \
    --dmnd_db eggnog_proteins.dmnd \
    --data_dir ./ \
    --cpu ${task.cpus} \
    -o ${specimen}.egm &&
    pigz -p ${task.cpus} ${specimen}.egm.emapper.annotations
    """
}

/*
// retrieve the eggnogmapper db
process eggnogMapperDownloadDB {
    container 'golob/eggnog-mapper:1.0.3__bcw.0.3.1B'
    label = 'multithread'
    //errorStrategy "retry"
    afterScript "rm -rf /tmp/enmdb/"

    output:
        file("eggnogDB.tgz") into eggnogdb_tgz
    
    """
    set -e
    
    mkdir -p /tmp/enmdb/ &&
    download_eggnog_data.py none -y --data_dir /tmp/enmdb/ &&
    tar -czvf /tmp/eggnogDB.tgz --directory /tmp/enmdb/ . &&
    ln -s /tmp/eggnogDB.tgz eggnogDB.tgz
    """
}

process eggnogMap {
    container 'golob/eggnog-mapper:1.0.3__bcw.0.3.1B'
    label 'mem_veryhigh'
    //errorStrategy "retry"

    input:
        set val(specimen), file(faa) from specimen_prokka_faa_ch
        file(eggnogdb_tgz)

    output:
        set val(specimen), file(" ${specimen}.egm.emapper.annotations.gz") into eggnogmap_ch

    """
    set -e 

    mkdir -p db/emdb/ &&
    tar -I pigz -C db/emdb/ -xf ${eggnogdb_tgz} ./eggnog.db ./eggnog_proteins.dmnd 
    emapper.py \
    -i ${faa} \
    -m diamond \
    --dmnd_db db/emdb/eggnog_proteins.dmnd \
    --data_dir db/emdb/ \
    --cpu ${task.cpus} \
    -o ${specimen}.egm &&
    pigz ${specimen}.egm.emapper.annotations
    """
        

}



// */