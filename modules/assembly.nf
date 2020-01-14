// Processes to perform de novo assembly and annotate those assembled sequences

// Assembly with metaspades
process metaspadesAssembly {
    container 'golob/spades:3.13.1__bcw.0.3.1'
    label 'mem_veryhigh'
    errorStrategy "retry"

    publishDir "${params.output_folder}/assembly/${specimen}", mode: 'copy'

    input:
        tuple specimen, file(R1), file(R2)
    
    output:
        tuple specimen, file("${specimen}.contigs.fasta"), file("${specimen}.scaffolds.fasta"), file("${specimen}.metaspades.log")
    
"""
set -e 

metaspades.py \
    --meta \
    --phred-offset ${params.phred_offset} \
    -1 ${R1} -2 ${R2} \
    -o . \
    -t ${task.cpus} -m ${task.memory.toMega() / 1024} | 
    tee -a ${specimen}.metaspades.log

ls -lahtr
mv contigs.fasta ${specimen}.contigs.fasta
mv scaffolds.fasta ${specimen}.scaffolds.fasta
"""
}


// Annotation with prokka

process prokkaAnnotate {
    container 'golob/prokka:1.1.14__bcw.0.3.1'
    label 'mem_veryhigh'
    errorStrategy "retry"
    publishDir "${params.output_folder}/prokka/${specimen}/", mode: 'copy'

    input:
        tuple val(specimen) , file(contigs), file(scaffolds), file(spades_log)
    
    output:
        tuple \
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
            file("prokka/${specimen}.txt.gz")

        tuple val(specimen), file("prokka/${specimen}.faa.gz")
    
"""
set -e 

prokka \
    --outdir prokka/ \
    --prefix ${specimen} \
    --cpus ${task.cpus} \
    --metagenome \
    --compliant \
    --centre ${params.centre} \
    ${scaffolds}
    
gzip prokka/${specimen}*
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
    set min_identity, file("mmseqs.${params.min_identity}.tsv.gz"), file("mmseqs.${params.min_identity}.rep.fasta.gz")
    

    afterScript "rm -r *"

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
