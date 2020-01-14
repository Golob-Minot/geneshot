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
    // errorStrategy "retry"
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