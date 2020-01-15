// Processes used for alignment of reads against gene databases

// Align each sample against the reference database of genes using DIAMOND
process makeDiamondDB {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    file fasta

    output:
    file "${fasta}.dmnd"

    """
    set -e
    diamond \
      makedb \
      --in ${fasta} \
      --db ${fasta}.dmnd \
      --threads ${task.cpus}
    """

}


// Align each sample against the reference database of genes using DIAMOND
process diamond {
    container "quay.io/fhcrc-microbiome/famli@sha256:25c34c73964f06653234dd7804c3cf5d9cf520bc063723e856dae8b16ba74b0c"
    label 'mem_veryhigh'
    // errorStrategy 'retry'
    
    input:
    tuple val(sample_name), file(R1), file(R2)
    file refdb
    
    output:
    tuple sample_name, file("${sample_name}.aln.gz")

    """
    set -e

    cat ${R1} ${R2} > query.fastq.gz

    diamond \
      blastx \
      --query query.fastq.gz \
      --out ${sample_name}.aln.gz \
      --threads ${task.cpus} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --min-score ${params.dmnd_min_score} \
      --query-cover ${params.dmnd_min_coverage} \
      --id ${params.dmnd_min_identity} \
      --top ${params.dmnd_top_pct} \
      --block-size ${task.memory.toMega() / (1024 * 6)} \
      --query-gencode ${params.gencode} \
      --compress 1 \
      --unal 0
    """

}


// Filter the alignments with the FAMLI algorithm
process famli {
    container "quay.io/fhcrc-microbiome/famli@sha256:241a7db60cb735abd59f4829e8ddda0451622b6eb2321f176fd9d76297d8c9e7"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    tuple sample_name, file(input_aln)
    
    output:
    tuple sample_name, "${sample_name}.json.gz"

    """
    set -e
    famli \
      filter \
      --input ${input_aln} \
      --output ${sample_name}.json \
      --threads ${task.cpus} \
      --batchsize 5000000
    gzip ${sample_name}.json
    """

}