#!/usr/bin/env python3

import os
import pandas as pd
import sys

# Set the input path
results_hdf = sys.argv[1]
print("Results HDF5: %s" % (results_hdf))

# Set the query string
query_str = sys.argv[2]
print("Query string: %s" % (query_str))

# Set the output path
output_fp = sys.argv[3]
print("Output filepath: %s" % (output_fp))

# Make sure that the file is present in the working folder
assert os.path.exists(results_hdf)

# Read in the table with the names of each eggNOG annotation
eggnog_desc = pd.read_hdf(
    results_hdf,
    "/ref/eggnog_desc"
).set_index(
    "eggnog_desc_ix"
)["eggnog_desc"]

# Filter the complete set of annotations to just those which
# contain the query string of interest
eggnog_desc = eggnog_desc.loc[eggnog_desc.fillna('').apply(lambda s: query_str in s)]

assert eggnog_desc.shape[0] > 0, "No annotations found"

print("There are %s annotations containing '%s':" % (str(eggnog_desc.shape[0]), query_str))
print('\n'.join(eggnog_desc.sort_values().tolist()))

# Set up a function to filter a DataFrame by the query string
def filter_df(df):

    # Add the annotation (if any)
    df = df.assign(
        eggNOG_desc=df['eggNOG_desc_ix'].fillna(-1).apply(int).apply(eggnog_desc.get)
    )

    # Filter out any genes which are missing annotations
    return df.reindex(index=df['eggNOG_desc'].dropna().index.values)

# Read in the table in chunks and filter as we go
df = pd.concat([
    filter_df(chunk_df)
    for chunk_df in pd.read_hdf(
        results_hdf, 
        '/annot/gene/all', 
        iterator=True
    )
])
print("Number of genes containing the query '%s': %d" % (query_str, df.shape[0]))

# If there are any genes matching this string
if df.shape[0] > 0:

    # Write out the smaller table
    df.to_csv(output_fp, index=None)

    print("Done")

else:

    print("NO GENES FOUND MATCHING THE QUERY: %s" % query_str)