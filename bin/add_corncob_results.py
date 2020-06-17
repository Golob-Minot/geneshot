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


# Read in the table of corncob results
corncob_wide = read_corncob_results()

# Read in the gene annotations
gene_annot = pd.read_hdf(hdf_fp, "/annot/gene/all")

# Open a connection to the HDF5
with pd.HDFStore(hdf_fp, "a") as store:

    # Write corncob results to HDF5
    corncob_wide.to_hdf(store, "/stats/cag/corncob")
