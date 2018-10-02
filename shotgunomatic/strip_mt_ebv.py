#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import sys

###
#   For shotgun microbiome experiments, we wish to strip out human
#   genomic DNA (for various reasons, including privacy concerns).
#   A typical pipeline will align a human reference genome to remove any
#   reads from human genomic DNA sequences.
#
#   The human reference genomes include not just human chromsomal DNA
#   but also human mitochondrial DNA and the EBV genome.
#   ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/
#
#   For microbiome studies, this isn't desirable.
#   Herpesvirii (including EBV, CMV, HHV6) are an important part
#   of the microbiome, as are bacterial genomes (including crucially the
#   16S rRNA gene) that can share components with huMT-DNA.
#
#   Thus, I have created this utility to take a human reference genome FASTA
#   file and remove the MT-DNA and EBV DNA (as identified by ID's below).
###

#  Which sequence IDs should be removed from this reference genome?
STRIPPED_IDS = {
    'chrM',  # Mitochondrial DNA. (Want to remove to ease 16S extraction)
    'chrEBV'  # EBV sequence
}


def main():
    """Entrypoint for main script."""
    arg_parser = argparse.ArgumentParser(description="""
    Remove mitochondrial and EBV DNA from the reference human genome.
    For gzipped sequences, consider using pipes: 
    gunzip -c ref.fna.gz | strip_mt_ebv.py | gzip > ref.nomtebv.fna.gz
    """)

    arg_parser.add_argument(
        '--infile', '-I',
        type=argparse.FileType('r'),
        default=sys.stdin
    )
    arg_parser.add_argument(
        '--outfile', '-O',
        type=argparse.FileType('w'),
        default=sys.stdout
    )
    args = arg_parser.parse_args()

    for sr in SeqIO.parse(args.infile, 'fasta'):
        if sr.id not in STRIPPED_IDS:
            SeqIO.write(sr, args.outfile, 'fasta')

if __name__ == "__main__":
    main()
