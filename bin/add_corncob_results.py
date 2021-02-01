#!/usr/bin/env python3

from functools import lru_cache
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import sys
import pickle
pickle.HIGHEST_PROTOCOL = 4

# Set the path to the HDF
hdf_fp = sys.argv[1]
corncob_csv = sys.argv[2]
fdr_method = sys.argv[3]
group_name = sys.argv[4]

def read_corncob_results():
    # Read in the corncob results
    corncob_df = pd.read_csv(corncob_csv)

    print(
        "Read in corncob results for %d %s groups" % 
        (corncob_df[group_name].unique().shape[0], group_name)
    )

    # Make a wide version of the table with just the mu. values
    corncob_wide = corncob_df.loc[
        corncob_df["parameter"].apply(
            lambda s: s.startswith("mu.")
        )
    ].pivot_table(
        index = [group_name, "parameter"],
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

    # Add the wald metric
    corncob_wide = corncob_wide.assign(
        wald=corncob_wide["estimate"] / corncob_wide["std_error"]
    )

    return corncob_wide


# Read in the table of corncob results
corncob_wide = read_corncob_results()

# Open a connection to the HDF5
with pd.HDFStore(hdf_fp, "a") as store:

    # Write corncob results to HDF5
    corncob_wide.to_hdf(store, "/stats/%s/corncob" % (group_name.lower()))
