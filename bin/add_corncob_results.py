#!/usr/bin/env python3

from functools import lru_cache
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import sys

# Set the path to the HDF
hdf_fp = sys.argv[1]
corncob_csv = sys.argv[2]
fdr_method = sys.argv[3]

def read_corncob_results(fdr_method=fdr_method, corncob_csv=corncob_csv):
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


# Taxonomy used in this analysis
class Taxonomy:
    def __init__(self, hdf_fp):

        # Read in the taxonomy table contained in the HDF store
        self.taxonomy_df = pd.read_hdf(
            hdf_fp, 
            "/ref/taxonomy"
        ).set_index(
            "tax_id"
        )

        # Make sure that the parent is encoded as a string
        #  (missing values otherwise would coerce it to a float)
        self.taxonomy_df["parent"] = self.taxonomy_df[
            "parent"
        ].fillna(
            0
        ).apply(
            float
        ).apply(
            int
        ).apply(
            str
        )

    def parent(self, tax_id):
        return self.taxonomy_df["parent"].get(tax_id)

    def ancestors(self, tax_id, a=[]):

        a = a + [tax_id]

        if self.parent(tax_id) is None:
            return a
        elif self.parent(tax_id) == tax_id:
            return a
        else:
            return self.ancestors(
                self.parent(tax_id),
                a=a.copy()
            )

    def name(self, tax_id):
        return self.taxonomy_df["name"].get(tax_id)

    def rank(self, tax_id):
        return self.taxonomy_df["rank"].get(tax_id)


def annotate_taxa(gene_annot, hdf_fp, rank_list=["family", "genus", "species"]):
    """Function to add species, genus, etc. labels to the gene annotation."""

    # Read in the taxonomy
    tax = Taxonomy(hdf_fp)

    # Add an LRU cache to the `ancestors` function
    @lru_cache(maxsize=None)
    def ancestors(tax_id):
        return tax.ancestors(str(tax_id))

    # Add an LRU cache to the `anc_at_rank` function
    @lru_cache(maxsize=None)
    def anc_at_rank(tax_id, rank="species"):
        if tax_id == 0:
            return
        for t in ancestors(str(tax_id)):
            if tax.rank(t) == rank:
                return tax.name(t)

    # Add the rank-specific taxon names
    for rank in rank_list:
        print("Adding {} names for genes".format(rank))

        gene_annot = gene_annot.assign(
            new_label=gene_annot["tax_id"].apply(lambda t: anc_at_rank(t, rank=rank))
        ).rename(
            columns={"new_label": rank}
        )

    print("Finished adding taxonomic labels")

    return gene_annot


def calc_enrichment(corncob_wide, gene_annot, col_name):
    """Function to calculate the enrichment of each label with each parameter from the geneshot output."""

    # Calculate the enrichment for each parameter individually
    results = pd.concat([
        calc_enrichment_parameter(
            parameter_name, 
            corncob_parameter_df.set_index(
                "CAG"
            ).reindex(
                columns=["p_value", "estimate"]
            ),
            gene_annot, 
            col_name
        )
        for parameter_name, corncob_parameter_df in corncob_wide.groupby("parameter")
    ])

    # Add the FDR corrected p-values
    results = results.assign(
        qvalue=multipletests(results["pvalue"], 0.2, "fdr_bh")[1]
    )

    return results


def calc_enrichment_parameter(parameter_name, cag_pvalues, gene_annot, col_name):

    print("Calculating enrichment for parameter {} by {}".format(
        parameter_name, col_name
    ))

    # Make sure that we have access to the indicated annotations
    assert col_name in gene_annot.columns.values, \
        "{} not a valid column name".format(col_name)

    # Subset to those genes belonging to CAGs which also have a p-value of any sort
    gene_annot_parameter = gene_annot.loc[
        gene_annot["CAG"].isin(cag_pvalues.index.values)
    ]
    print("There are {:,} genes from {:,} CAGs with p-values for {}".format(
        gene_annot_parameter.shape[0], gene_annot_parameter["CAG"].unique().shape[0], parameter_name
    ))

    assert gene_annot_parameter.shape[0] > 0

    # Treat the CAGs separately which have positive and negative estimated coefficients
    pos_coef_pvalues = cag_pvalues.query("estimate > 0")["p_value"]
    neg_coef_pvalues = cag_pvalues.query("estimate < 0")["p_value"]

    # Keep track of all of the results
    results = []

    # For each annotation, pick which CAGs have it and which don't
    for annotation_label, annotated_genes in gene_annot.groupby(col_name):
        # If there is just a single CAG with this annotation, skip it
        if annotated_genes["CAG"].unique().shape[0] == 1:
            continue

        # Consider separately the CAGs with estimated coefficients > 0 or < 0
        for subset_label, cag_pvalues_subset in [
            ("positive", pos_coef_pvalues),
            ("negative", neg_coef_pvalues),
        ]:

            # Get the p-values for CAGs with and without this annotation
            pvalues_with_annot = cag_pvalues_subset.reindex(
                index=annotated_genes["CAG"].unique()).dropna()
            pvalues_lacking_annot = cag_pvalues_subset.drop(
                index=annotated_genes["CAG"].unique(),
                errors="ignore"
            ).dropna()

            if pvalues_with_annot.shape[0] == 0 or pvalues_lacking_annot.shape[0] == 0:
                continue

            statistic, pvalue = stats.mannwhitneyu(
                pvalues_with_annot.values,
                pvalues_lacking_annot.values,
                alternative="less"
            )

            results.append({
                "label": annotation_label,
                "category": col_name,
                "n_with_label": pvalues_with_annot.shape[0],
                "n_lacking_label": pvalues_lacking_annot.shape[0],
                "statistic": statistic,
                "pvalue": pvalue,
                "est_coef": subset_label,
                "parameter": parameter_name
            })

    print("Processed {:,} {} labels for {}".format(len(results), col_name, parameter_name))

    return pd.DataFrame(results)


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
