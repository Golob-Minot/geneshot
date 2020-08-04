container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.0.3"

process linclust {
    tag "Cluster genes with similar sequences"
    container "quay.io/fhcrc-microbiome/integrate-metagenomic-assemblies:v0.5"
    label 'mem_medium'
    errorStrategy 'retry'
    
    input:
    file "input.genes.*.fasta.gz"
    
    output:
    file "output.genes.fasta.gz"
    
"""
#!/bin/bash

set -e

# Combine input files
echo "Combining input files"
cat input.genes.* > input.genes.fasta.gz

# Make the MMSeqs2 database
echo "Running linclust"
mmseqs createdb input.genes.fasta.gz db

# Cluster the protein sequences
mmseqs linclust db cluster_db ./ \
    --min-seq-id ${params.min_identity / 100} \
    --max-seqs 100000 \
    -c ${params.min_coverage / 100}

# Get the representative sequences
mmseqs result2repseq db cluster_db genes
mmseqs result2flat db db genes output.genes.fasta --use-fasta-header
echo "Compressing"
gzip output.genes.fasta

echo "Done"
"""
}


process diamondDedup {
    tag "Deduplicate sequences by alignment with DIAMOND"
    container "quay.io/fhcrc-microbiome/famli:v1.5"
    label 'mem_veryhigh'
    errorStrategy 'retry'
    
    input:
    file "input.genes.fasta.gz"
    
    output:
    file "output.genes.fasta.gz"
    
"""
#!/bin/bash

set -e

# Make a DIAMOND database
diamond \
    makedb \
    --in input.genes.fasta.gz \
    --db input.genes.dmnd \
    --threads ${task.cpus}

# Align the genes against themselves and filter any which align
diamond \
    blastp \
    --query input.genes.fasta.gz \
    --threads ${task.cpus} \
    --db input.genes.dmnd \
    --out input.genes.aln \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
    --query-cover 90 \
    --id ${params.min_identity} \
    --top 0 \
    --block-size ${task.memory.toMega() / (1024 * 6)} \
    --unal 0

python << END

from Bio.SeqIO.FastaIO import SimpleFastaParser
import gzip
import sys

# Keep track of the genes which have been filtered out
duplicate_genes = set([])

# Iterate over the alignment
print("Reading alignments")
ix = 0
for line in open("input.genes.aln", "r"):

    ix += 1

    if ix % 100000 == 0:
        print("Read %d lines of alignment - found %d duplicated genes" % (ix, len(duplicate_genes)))

    qname, sname, _, _, _, _, _, _, _, _, _, _, qlen, slen = line.rstrip("\\n").split("\\t")

    # Skip self alignments
    if qname == sname:
        continue
    # If we have excluded either of the genes before, skip the line
    if qname in duplicate_genes or sname in duplicate_genes:
        continue
    # For non-self alignments, remove the smaller of the two
    if int(slen) < int(qlen):
        duplicate_genes.add(sname)
    else:
        duplicate_genes.add(qname)

assert ix > 0, "Didn't read any alignments"
print("Read %d lines of alignment - found %d duplicated genes" % (ix, len(duplicate_genes)))
print("Done reading alignments")

# Now let's make the filtered FASTA with all duplicated genes removed
n_found = 0
n = 0
with gzip.open("input.genes.fasta.gz", "rt") as handle_in:
    with gzip.open("output.genes.fasta.gz", "wt") as handle_out:
        for header, seq in SimpleFastaParser(handle_in):

            header = header.split(" ", 1)[0]

            n += 1
            if header in duplicate_genes:
                n_found += 1
            else:
                handle_out.write(">%s\\n%s\\n" % (header, seq))

# Make sure that we encountered all of the duplicated genes
print("Read in %d sequences, filtered out %d, wrote out the rest" % (n, n_found))
assert n_found == len(duplicate_genes), "%d != %d" % (n_found, len(duplicate_genes))

END

echo "Done"
"""
}
