#!/usr/bin/env python3

import pandas as pd
from statsmodels.stats.multitest import multipletests

# Set the path to the HDF
hdf_fp = "${results_hdf}"


def read_corncob_results(fdr_method="${params.fdr_method}", corncob_csv= "${corncob_csv}"):
    # Read in the corncob results
    corncob_df = pd.read_csv(corncob_csv)

    print(
        "Read in corncob results for %d CAGs" % 
        corncob_df["CAG"].unique().shape[0]
    )

    # Make a wide version of the table with just the mu. values
    corncob_wide = corncob_df.loc[
        corncob_df["parameter"].apply(
            lambda s: s.startswith("mu.")
        )
    ].pivot_table(
        index = ["CAG", "parameter"],
        columns = "type",
        values = "value"
    ).reset_index(
    ).apply(
        lambda v: v.apply(lambda s: s.replace("mu.", "")) if v.name == "parameter" else v
    )

    # Adding the q-value is conditional on p-values being present
    if "p_value" in corncob_wide.columns.values:

        # Add the q-value (FDR-BH)
        corncob_wide = corncob_wide.assign(
            q_value = multipletests(corncob_wide.p_value.fillna(1), 0.2, fdr_method)[1]
        )

    return corncob_wide


def annotate_taxa(gene_annot, hdf_fp):
    """Function to add species, genus, etc. labels to the gene annotation."""
    return gene_annot


def calc_enrichment(corncob_wide, gene_annot, rank):
    """Function to calculate the enrichment of each label with each parameter from the geneshot output."""
    return


# Read in the table of corncob results
corncob_wide = read_corncob_results()

# Read in the gene annotations
gene_annot = pd.read_hdf(hdf_fp, "/annot/gene/all")

# If we are able to calculate any label enrichment scores, save them in this dict
enrichment_dict = dict()

# Check if we have taxonomic assignments annotating the gene catalog
if "tax_id" in gene_annot.columns.values:
    gene_annot = annotate_taxa(gene_annot, hdf_fp)

    for rank in ["species", "genus", "family"]:
        enrichment_dict[rank] = calc_enrichment(corncob_wide, gene_annot, rank)

# Check if we have eggNOG annotations for each gene
if "eggNOG_desc" in gene_annot.columns.values:
    enrichment_dict["eggNOG_desc"] = calc_enrichment(corncob_wide, gene_annot, "eggNOG_desc")

# Open a connection to the HDF5
with pd.HDFStore(hdf_fp, "a") as store:

    # Write corncob results to HDF5
    corncob_wide.to_hdf(store, "/stats/cag/corncob")

    # Write any enrichment results to the HDF
    if len(enrichment_dict) > 0:
        for enrichment_label, enrichment_df in enrichment_dict.items():
            if enrichment_df is not None:
                enrichment_df.to_hdf(store, "/stats/enrichment/{}".format(enrichment_label))
