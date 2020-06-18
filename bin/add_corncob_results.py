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


def write_corncob_by_annot(corncob_wide, gene_annot, col_name_list, fp_out_template):
    
    for col_name in col_name_list:
        assert col_name in gene_annot.columns.values

    # Write out chunks of data
    batch_size = 0
    ix = 0
    batch = []

    for df in [
        corncob_wide.loc[
            corncob_wide["CAG"].isin(genes_with_label["CAG"].unique())
        ].assign(
            label = label
        ).assign(
            annotation = col_name
        )
        for col_name in col_name_list
        for label, genes_with_label in gene_annot.groupby(col_name)
    ]:
        batch_size += df.shape[0]
        batch.append(df)

        if batch_size >= 100:
            pd.concat(batch).to_csv(
                fp_out_template.format(ix),
                index=None
            )
            ix += 1
            batch_size = 0
            batch = []


# Read in the table of corncob results
corncob_wide = read_corncob_results()

# Read in the gene annotations
gene_annot = pd.read_hdf(hdf_fp, "/annot/gene/all")

# Get the list of columns to use
columns_to_use = []

# Check if we have taxonomic assignments annotating the gene catalog
if "tax_id" in gene_annot.columns.values:

    # Add the names for the species, genus, and family
    # which will become new columns in `gene_annot` 
    gene_annot = annotate_taxa(gene_annot, hdf_fp)

    for rank in ["species", "genus", "family"]:

        columns_to_use.append(rank)

# Check if we have eggNOG annotations for each gene
if "eggNOG_desc" in gene_annot.columns.values:

    columns_to_use.append("eggNOG_desc")

# Write out a CSV with corncob results for each
# CAG which are grouped according to which CAGS
# have a given annotation
# This will be used to run betta in a separate step
write_corncob_by_annot(
    corncob_wide,
    gene_annot,
    columns_to_use,
    "corncob.shard.{}.csv.gz"
)

# Open a connection to the HDF5
with pd.HDFStore(hdf_fp, "a") as store:

    # Write corncob results to HDF5
    corncob_wide.to_hdf(store, "/stats/cag/corncob")
