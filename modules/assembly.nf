// Processes to perform de novo assembly and annotate those assembled sequences

workflow assembly_wf {
    get:
        combined_reads_ch

    main:

    // Assemble samples with metaSPAdes
    metaspadesAssembly(
        combined_reads_ch
    )

    // Annotate those contigs with Prokka
    prodigalAnnotate(
        metaspadesAssembly.out
    )

    // Combine the gene sequences across all samples
    combineCDS(
        prodigalAnnotate.out[
            0
        ].map {
            it -> it[1]
        }.collect()
    )

    // Combine genes by amino acid identity
    clusterCDS(
        combineCDS.out
    )

    // Determine whether or not to run the eggNOG annotation based
    // on --noannot and --eggnog_db / --eggnog_dmnd
    run_eggnog = false
    if ( params.noannot == false ) {
        if ( params.eggnog_db && params.eggnog_dmnd ) {
            if ( !file(params.eggnog_db).isEmpty() && !file(params.eggnog_dmnd).isEmpty() ){
                run_eggnog = true
            }
        }
    }

    // Annotate the clustered genes with eggNOG
    if ( run_eggnog ){
        eggnog_annotation(
            clusterCDS.out[0],
            file(params.eggnog_db),
            file(params.eggnog_dmnd)
        )
    }

    // Determine whether or not to run the taxnomic annotation based
    // on --noannot and --taxonomic_dmnd
    run_tax = false
    if ( params.noannot == false ) {
        if ( params.taxonomic_dmnd ) {
            if ( !file(params.taxonomic_dmnd).isEmpty() ){
                run_tax = true
            }
        }
    }

    // Annotate the clustered genes with DIAMOND for taxonomic ID
    if ( run_tax ) {
        taxonomic_annotation(
            clusterCDS.out[0],
            path(params.taxonomic_dmnd)
        )
    }

    emit:
        gene_fasta = clusterCDS.out[0]

}

// Assembly with metaspades
process metaspadesAssembly {
    container 'golob/spades:3.13.1__bcw.0.3.1'
    label 'mem_veryhigh'
    errorStrategy "retry"

    publishDir "${params.output_folder}/assembly/${specimen}", mode: "copy"

    input:
        tuple specimen, file(R1), file(R2)
    
    output:
        tuple specimen, file("${specimen}.contigs.fasta.gz"), file("${specimen}.scaffolds.fasta.gz"), file("${specimen}.metaspades.log")
    
"""
set -e 

metaspades.py \
    --meta \
    --phred-offset ${params.phred_offset} \
    -1 ${R1} -2 ${R2} \
    -o . \
    -t ${task.cpus} -m ${task.memory.toMega() / 1024} | 
    tee -a ${specimen}.metaspades.log

gzip -c contigs.fasta > ${specimen}.contigs.fasta.gz
gzip -c scaffolds.fasta > ${specimen}.scaffolds.fasta.gz
"""
}

// Annotation of coding sequences with prodigal

process prodigalAnnotate {
    container 'quay.io/biocontainers/prodigal:2.6.3--h516909a_2'
    label 'io_limited'
    errorStrategy "retry"
    publishDir "${params.output_folder}/prodigal/${specimen}/", mode: "copy"

    input:
        tuple val(specimen), file(contigs), file(scaffolds), file(spades_log)
    
    output:
        tuple val(specimen), file("${specimen}.faa.gz")
        tuple val(specimen), file("${specimen}.gff.gz")
    
"""
set -e 

gunzip -c ${contigs} > ${specimen}.contigs.fasta

prodigal \
    -a ${specimen}.faa \
    -i  ${specimen}.contigs.fasta \
    -f gff \
    -o ${specimen}.gff \
    -p meta

gzip ${specimen}.gff
gzip ${specimen}.faa

"""
}

process combineCDS {
    container "ubuntu:16.04"
    label 'io_limited'
    errorStrategy 'retry'
    maxRetries 10
    
    input:
    file "assembly.*.fasta.gz"
    
    output:
    file "all_CDS.fasta.gz"

    """
#!/bin/bash
set -e
cat *fasta.gz > all_CDS.fasta.gz
gzip -t all_CDS.fasta.gz
    """
}


process clusterCDS {
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    file all_cds
    
    output:
    file "mmseqs.${params.min_identity}.rep.fasta.gz"
    file "mmseqs.${params.min_identity}.tsv.gz"
    
    """
#!/bin/bash
set -e
# Make the MMSeqs2 database
mmseqs createdb ${all_cds} db
# Cluster the protein sequences
mmseqs linclust db mmseqs.${params.min_identity}.cluster ./ \
    --min-seq-id ${params.min_identity / 100} \
    --max-seqs 100000 \
    -c ${params.min_coverage / 100}
# Make TSV output for clustering
mmseqs createtsv db db mmseqs.${params.min_identity}.cluster mmseqs.${params.min_identity}.tsv
# Get the representative sequences
mmseqs result2repseq db mmseqs.${params.min_identity}.cluster mmseqs.${params.min_identity}.rep
mmseqs result2flat db db mmseqs.${params.min_identity}.rep mmseqs.${params.min_identity}.rep.fasta --use-fasta-header
gzip mmseqs.${params.min_identity}.tsv
gzip mmseqs.${params.min_identity}.rep.fasta
    """
}


process taxonomic_annotation {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'

    input:
    file query
    file diamond_tax_db
    
    output:
    file "${query}.tax.aln.gz"

    
"""
set -e

diamond \
    blastp \
    --db ${diamond_tax_db} \
    --query ${query} \
    --out ${query}.tax.aln.gz \
    --outfmt 102 \
    --id ${params.min_identity} \
    --top ${100 - params.min_identity} \
    --block-size ${task.memory.toMega() / (1024 * 6)} \
    --threads ${task.cpus} \
    --compress 1

rm ${diamond_tax_db}
"""

}


process eggnog_annotation {
    container "quay.io/biocontainers/eggnog-mapper:2.0.1--py_1"
    label 'mem_veryhigh'
    
    input:
    path query
    path eggnog_db
    path eggnog_dmnd

    output:
    path "${query}.emapper.annotations.gz"

    
    """
set -e

mkdir data
mkdir TEMP
mkdir SCRATCH

mv ${eggnog_db} data/eggnog.db
mv ${eggnog_dmnd} data/eggnog_proteins.dmnd

emapper.py \
    -i ${query} \
    --output ${query} \
    -m "diamond" \
    --cpu ${task.cpus} \
    --data_dir data/ \
    --scratch_dir SCRATCH/ \
    --temp_dir TEMP/ \

gzip ${query}.emapper.annotations
    
    """

}

