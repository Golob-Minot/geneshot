#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from scipy.spatial import distance
from matplotlib.backends.backend_pdf import PdfPages
from  matplotlib import pyplot as plt
import seaborn as sns


def summarize_cag_coabundance(
    results_hdf=None,
    details_hdf=None,
    output=None,
    metric="cosine",
    n=1000
):

    assert os.path.exists(results_hdf), "Please specify --results-hdf"
    assert os.path.exists(details_hdf), "Please specify --details-hdf"
    assert output is not None, "Please specify --output"

    # Read in the table of CAG assignments
    cag_df = pd.read_hdf(
        results_hdf,
        "/annot/gene/all"
    ).reindex(
        columns=["gene", "CAG"]
    )

    # Get the unique set of CAGs
    cag_list = cag_df["CAG"].drop_duplicates()
    
    # Pick a random set of CAGs
    cags_to_sample = cag_list.sample(
        min(n, cag_list.shape[0])
    )

    # Select two genes from each CAG
    cag_df = cag_df.loc[
        cag_df["CAG"].isin(cags_to_sample).values
    ].groupby("CAG").head(2)

    # Set up a dict with the CAG for each gene
    cag_dict = cag_df.set_index("gene")["CAG"]

    # Make an abundance matrix for those genes
    abund_df = make_abund_df(
        details_hdf,
        cag_df["gene"].tolist()
    )

    # Calculate the pairwise distance for all genes
    dists = pd.DataFrame(
        distance.squareform(
            distance.pdist(
                abund_df,
                metric=metric
            )
        ),
        index=abund_df.index.values,
        columns=abund_df.index.values
    )

    # Count up the distances for each pair
    output_df = pd.DataFrame([
        {
            "same_cag": cag_dict[gene_a] == cag_dict[gene_b],
            "distance": dists.loc[gene_a, gene_b]
        }
        for gene_a in abund_df.index.values
        for gene_b in abund_df.index.values
        if gene_a < gene_b
    ])

    with PdfPages(f"{output}.pdf") as pdf:

        g = sns.catplot(
            data=output_df,
            x="same_cag",
            y="distance",
            kind="box"
        )

        pdf.savefig()

        plt.close()

# Get the list of specimens which have assembly information
def get_specimen_list(details_hdf, prefix="/abund/gene/long/"):
    with pd.HDFStore(details_hdf, 'r') as store:
        return [
            p[len(prefix):]
            for p in store
            if p.startswith(prefix)
        ]


def read_abund(store, specimen, gene_id_list):

    # Read in the abundances
    df = pd.read_hdf(store, f"/abund/gene/long/{specimen}")

    # Get the sum of abundances
    abund_sum = df["depth"].sum()

    # Return the abundances for just the set of specified genes
    return df.set_index(
        "id"
    ).reindex(
        gene_id_list
    )[
        "depth"
    ].fillna(
        0
    ) / abund_sum

def make_abund_df(details_hdf, gene_id_list):

    # Get the list of specimens
    specimen_list = get_specimen_list(details_hdf)

    with pd.HDFStore(details_hdf, 'r') as store:
        
        abund_df = pd.DataFrame({
            specimen: read_abund(
                store,
                specimen,
                gene_id_list
            )
            for specimen in specimen_list
        })


    return abund_df

if __name__ == "__main__":

    ###################
    # PARSE ARGUMENTS #
    ###################

    # Create the parser
    parser = argparse.ArgumentParser(
        description='Summarize CAG co-abundance values'
    )

    # Add the arguments
    parser.add_argument(
        '--results-hdf',
        type=str,
        help='Geneshot output - results.hdf5'
    )
    # Add the arguments
    parser.add_argument(
        '--details-hdf',
        type=str,
        help='Geneshot output - results.hdf5'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Prefix for output files (*.csv and *.pdf)'
    )
    parser.add_argument(
        '--metric',
        type=str,
        default='cosine',
        help='Distance metric to be used for clustering'
    )
    parser.add_argument(
        '--n',
        type=int,
        default=1000,
        help='Number of CAGs to analyze (1 pair of genes each)'
    )

    # Parse the arguments
    args = parser.parse_args()

    summarize_cag_coabundance(
        **args.__dict__
    )