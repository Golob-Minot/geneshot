#!/usr/bin/env python3

import argparse
from collections import defaultdict
import gzip
import os
import pandas as pd


def split_fasta_CAGs(results_hdf=None, fasta=None, output_folder=None, output_prefix=None):
    """Split up a FASTA by CAG assignment."""

    assert results_hdf is not None, "Please provide --results-hdf"
    assert os.path.exists(results_hdf), f"Not found: {results_hdf}"
    assert fasta is not None, "Please provide --fasta"
    assert os.path.exists(fasta), f"Not found: {fasta}"
    assert output_prefix is not None, "Please provide --output-prefix"
    assert output_folder is not None, "Please provide --output-folder"
    if not os.path.exists(output_folder):
        print(f"Making output folder {output_folder}")
        os.mkdir(output_folder)

    # Open a connection to the HDF store
    with pd.HDFStore(results_hdf, "r") as store:

        # Read in the table of gene annotations
        annot_df = pd.read_hdf(store, "/annot/gene/all")

        # If we have eggNOG annotation data available
        if "/ref/eggnog_desc" in store:
            eggnog_name_dict = pd.read_hdf(
                results_hdf,
                "/ref/eggnog_desc"
            ).set_index(
                "eggnog_desc_ix"
            )[
                "eggnog_desc"
            ]

        else:
            eggnog_name_dict = dict()

    # Format a name for each gene
    gene_names = annot_df.set_index(
        "gene"
    ).reindex(
        columns=["tax_name", "eggNOG_desc_ix"]
    ).apply(
        lambda r: format_gene_name(r, eggnog_name_dict),
        axis=1
    )

    # Format the CAG assignment for each gene
    cag_assignment = annot_df.reindex(
        columns=["gene", "CAG"]
    ).dropna(
    ).set_index(
        "gene"
    )[
        "CAG"
    ].apply(int).to_dict()

    # Set up a dict with a list of tuples for each CAG
    cags = defaultdict(list)

    # Open the FASTA file
    with gzip.open(fasta, "rt") as handle:

        for header, seq in parse_fasta(handle):

            # Get the CAG for this sequence
            cag_id = cag_assignment.get(header)

            # If there is a CAG assigned
            if cag_id is not None:

                # Add this gene to the list
                cags[cag_id].append(
                    (gene_names.get(header), seq)
                )

    # Iterate over each CAG
    for cag_id, gene_list in cags.items():

        output_fp = os.path.join(
            output_folder,
            f"{output_prefix}{cag_id}.fasta.gz"
        )

        print(f"Writing out to {output_fp}")

        # Open a path to write out
        with gzip.open(output_fp, 'wt') as handle_out:

            handle_out.write("\n".join([
                f">{header}\n{seq}"
                for header, seq in gene_list
            ]))

def parse_fasta(handle):
    header = None
    seq = ""

    for line in handle:

        line = line.rstrip('\n')

        if line[0] == ">":

            if header is not None:
                yield header, seq

            header = line[1:].split(" ")[0]
            seq = ""

        else:
            seq = f"{seq}{line}"

    if header is not None:
        yield header, seq

def format_gene_name(r, eggnog_name_dict):

    # Start with the ID of the gene
    gene_name = r.name

    # If there is a taxonomic assignment
    if not pd.isnull(r["tax_name"]):

        # Add that taxonomic name to the gene
        gene_name = f"{gene_name} ({r['tax_name']})"

    # If there is a functional annotation
    if not pd.isnull(r.eggNOG_desc_ix):

        # Get the name for that function
        func_name = eggnog_name_dict.get(int(r.eggNOG_desc_ix))

        if not pd.isnull(func_name):

            gene_name = f"{gene_name} {func_name}"

    return gene_name

# Entrypoint
if __name__ == "__main__":
        
    ###################
    # PARSE ARGUMENTS #
    ###################

    # Create the parser
    parser = argparse.ArgumentParser(
        description='Split up a FASTA by CAG assignment'
    )

    # Add the arguments
    parser.add_argument(
        '--results-hdf',
        type=str,
        help='Geneshot results in HDF5 format'
    )
    # Add the arguments
    parser.add_argument(
        '--fasta',
        type=str,
        help='Gene sequences in FASTA format'
    )
    parser.add_argument(
        '--output-folder',
        type=str,
        default='./',
        help='Folder for output files'
    )
    parser.add_argument(
        '--output-prefix',
        type=str,
        default="CAG",
        help='Prefix for output files'
    )

    # Parse the arguments
    args = parser.parse_args()

    split_fasta_CAGs(
        results_hdf=args.results_hdf,
        fasta=args.fasta,
        output_folder=args.output_folder,
        output_prefix=args.output_prefix,
    )