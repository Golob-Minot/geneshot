#!/usr/bin/env python3

import argparse
import logging
import os
import pandas as pd

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [validate_geneshot_output.py] %(message)s'
)
rootLogger = logging.getLogger()
rootLogger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
rootLogger.addHandler(consoleHandler)


def validate_results_hdf(results_hdf, check_corncob = False):
    """Validate that the results HDF has all expected data."""

    for key_name in [
        "/manifest",
        "/summary/all",
        "/annot/gene/cag",
        "/annot/gene/all",
        "/annot/cag/all",
        "/abund/cag/wide",
        "/ordination/pca",
        "/ordination/tsne",
        "/distances/euclidean",
        "/distances/braycurtis",
        "/distances/jaccard",
        "/distances/aitchison",
    ]:
        # Try to read the table
        df = read_table_from_hdf(results_hdf, key_name)

        # Report errors with empty rows
        assert df.shape[0] > 0, "{} in {} has 0 rows".format(key_name, results_hdf)

    if check_corncob:

        # Read the results from running corncob
        df = read_table_from_hdf(results_hdf, "/stats/cag/corncob")

        # Report errors with empty rows
        assert df.shape[0] > 0, "{} in {} has 0 rows".format(key_name, results_hdf)
        

    # Some tables should be indexed -- let's check for that
    for key_name, where_str in [
        ("/abund/cag/wide", "CAG == 0"),
        ("/abund/cag/wide", "CAG == 1"),
        ("/annot/gene/all", "CAG == 0"),
        ("/annot/gene/all", "CAG == 1"),
    ]:
        df = read_table_from_hdf(results_hdf, key_name, where = where_str)
        msg = "Could not read index ({}) in {}".format(key_name, where_str)
        assert df.shape[0] > 0, msg


def validate_details_hdf(details_hdf, manifest_df, skip_assembly=False):
    """Validate that the details HDF has all expected data."""

    # Iterate over every sample from the manifest
    for sample_name in manifest_df["specimen"].values:
        # Iterate over every table which should be present
        for key_name in [
            "/abund/allele/assembly/",
            "/abund/gene/long/",
        ]:
            if key_name == "/abund/allele/assembly/" and skip_assembly:
                logging.info("Skipping assembly information (--skip-assembly was set)")
                continue

            # Try to read the table
            df = read_table_from_hdf(details_hdf, key_name + sample_name)

            # Report errors with empty rows
            msg = "{} in {} has 0 rows".format(key_name + sample_name, details_hdf)
            assert df.shape[0] > 0, msg


def read_table_from_hdf(hdf_fp, key_name, **kwargs):
    logging.info("Attempting to read {} from {} {}".format(
        key_name,
        hdf_fp,
        "({})".format(", ".join(["{}={}".format(k, v) for k, v in kwargs.items()])) if len(kwargs) > 0 else ""
    ))
    with pd.HDFStore(hdf_fp, "r") as store:
        assert key_name in store, "{} not found in {}".format(key_name, hdf_fp)
        return pd.read_hdf(store, key_name, **kwargs)


def validate_geneshot_output(results_hdf, details_hdf, skip_assembly = False, check_corncob = False):
    """Validate that the geneshot outputs have all expected data."""
    assert os.path.exists(results_hdf)
    assert os.path.exists(details_hdf)


    # Validate the results HDF
    validate_results_hdf(results_hdf, check_corncob=check_corncob)

    # Read the manifest from the results HDF
    manifest_df = read_table_from_hdf(results_hdf, "/manifest")

    # Validate the details HDF
    validate_details_hdf(details_hdf, manifest_df, skip_assembly=skip_assembly)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--results-hdf", 
        help = "Geneshot results output file (*.results.hdf5)",
        required = True
    )
    parser.add_argument(
        "--details-hdf", 
        help = "Geneshot details output file (*.details.hdf5)",
        required = True
    )
    parser.add_argument(
        "--skip-assembly", 
        help = "If specified, skip the detailed assembly information",
        action = "store_true"
    )
    parser.add_argument(
        "--check-corncob", 
        help = "If specified, check for corncob output",
        action = "store_true"
    )

    args = parser.parse_args()

    validate_geneshot_output(
        args.results_hdf,
        args.details_hdf,
        skip_assembly = args.skip_assembly,
        check_corncob = args.check_corncob
    )
