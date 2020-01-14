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
    errorStrategy 'retry'
    
    input:
    set val(sample_name), file(input_fastq)
    file refdb
    val min_id from 90
    val query_cover from 50
    val cpu from 32
    val top from 1
    val min_score from 20
    val blocks from 15
    val query_gencode from 11

    output:
    set sample_name, file("${sample_name}.aln.gz")

    afterScript "rm *"

    """
    set -e
    diamond \
      blastx \
      --query ${input_fastq} \
      --out ${sample_name}.aln.gz \
      --threads ${cpu} \
      --db ${refdb} \
      --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
      --min-score ${min_score} \
      --query-cover ${query_cover} \
      --id ${min_id} \
      --top ${top} \
      --block-size ${blocks} \
      --query-gencode ${query_gencode} \
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
    set sample_name, file(input_aln)
    val cpu from 16
    val batchsize from 50000000

    output:
    file "${sample_name}.json.gz"

    afterScript "rm *"

    """
    set -e
    famli \
      filter \
      --input ${input_aln} \
      --output ${sample_name}.json \
      --threads ${cpu} \
      --batchsize ${batchsize}
    gzip ${sample_name}.json
    """

}