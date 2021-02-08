#!/usr/bin/env python3

import argparse
from collections import defaultdict
from direct_redis import DirectRedis
from fastcluster import linkage
from functools import lru_cache
import logging
import math
import numpy as np
import os
import pandas as pd
from scipy.cluster.hierarchy import cophenet, optimal_leaf_ordering

#################
# DOCUMENTATION #
#################

# The purpose of this script is to index a set of geneshot results
# in a format which supports rapid retrieval for interactive visualization.
# The storage format we have selected is redis, a key-value store which
# supports Python objects. All data objects will be accessed under a single
# key, laid out as follows. Note that the <X> syntax indicates string
# interpolation.

# EXPERIMENT
# manifest
#   DataFrame with specimen, and additional user-provided metadata columns

# CAG METRICS
# cag_summary_metrics
#   DataFrame with cag, size, prevalence, mean_abundance, std_abundance, entropy

#  TAXONOMIC ANNOTATION
# taxonomy_df
#   DataFrame with tax_id, name, rank, parent
# taxonomy_name
#   Dict with <tax_id>:<tax_name>
# taxonomy_rank
#   Dict with <tax_id>:<tax_rank>
# cag_majority_taxon <phylum | class | order | family | genus | species>
#   Dict with <cag_id>:<name of majority taxon at the indicated rank>
# cag_tax_assignments <cag_id>
#   DataFrame with the number of genes assigned to each tax ID
# tax_cag_assignments <TAX_ID>
#   Dict with <cag_id>:<number of genes from CAG assigned to TAX_ID>

# FUNCTIONAL ANNOTATION
# func_index
#   Dict with <func_ix>:<long name of eggNOG function>
# cag_func_assignments <cag_id>
#   Dict with <func_ix>:<n_genes>
# func_cag_set <func_id>
#   Set of {all <cag_id> which contain at least one gene assigned}

# ASSOCIATIONS (CORNCOB)
# available_parameters
#   List with all <PARAM> values available
# cag_association <PARAM>
#   DataFrame with cag_id, estimate, p_value, std_error, q_value, wald
# tax_association <PARAM>
#   DataFrame with tax_id, estimate, p_value, std_error, q_value, wald
# func_association <PARAM>
#   DataFrame with func_id, estimate, p_value, std_error, q_value, wald

# ABUNDANCES
# cag_abundance_specimen <specimen>
#   Dict with <cag_id>:<relative abundance of CAG>
# specimen_abundance_cag <cag_id>
#   Dict with <specimen>:<relative abundance of CAG>
# specimen_abundance_func <func_id>
#   Dict with <specimen>:<relative abundance>
# func_abundance_specimen <specimen>
#   Dict with <func_id>:<relative abundance>
# specimen_abundance_tax <tax_id>
#   Dict with <specimen>:<relative abundance>
# tax_abundance_specimen <specimen>
#   Dict with <tax_id>:<relative abundance>

# ORDINATION
# specimen_ordination_pca <cag | taxa | func>
#   DataFrame with specimen, PC1 (X%), PC1 (X%), etc. 
# cag_taxonomic_layout
#   DataFrame with node, parent, terminal, x, y

# MEAN ABUNDANCES
# mean_abundance_cags
#   Series with <cag_id>:<mean abundance> in descending order
# mean_abundance_taxa
#   Series with <tax_id>:<mean abundance> in descending order
# mean_abundance_func
#   Series with <func_id>:<mean abundance> in descending order

##################
# SET UP LOGGING #
##################

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [buildRedis] %(message)s'
)
logger = logging.getLogger('buildRedis')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

###################
# PARSE ARGUMENTS #
###################

# Create the parser
parser = argparse.ArgumentParser(
    description='Format data for rapid visualization in a redis store'
)

# Add the arguments
parser.add_argument(
    '--results',
    type=str,
    help='Primary geneshot output (*.results.hdf5)'
)
parser.add_argument(
    '--details',
    type=str,
    help='Supplementary geneshot output (*.details.hdf5)'
)
parser.add_argument(
    '--host',
    type=str,
    default="localhost",
    help='Redis host'
)
parser.add_argument(
    '--port',
    type=int,
    default=6379,
    help='Redis port'
)

# Parse the arguments
args = parser.parse_args()


def buildRedis(results=None, details=None, host=None, port=None):
    """Reformat geneshot outputs to support rapid visualization."""

    logger.info(f"Making sure that {results} exists")
    assert os.path.exists(results)
    logger.info(f"Making sure that {details} exists")
    assert os.path.exists(details)

    # Connect to redis
    logger.info(f"Connecting to redis at {host}:{port}")
    r = DirectRedis(host=host, port=port)

    # Connect to the HDF input files
    logger.info("Opening HDF stores")
    with pd.HDFStore(results, 'r') as results_store, pd.HDFStore(details, 'r') as details_store:

        # Save basic CAG information
        save_cag_data(r, results_store, details_store)

        # Save experimental metadata
        save_expt_data(r, results_store, details_store)

        # Save taxonomic information
        save_tax_data(r, results_store, details_store)

        # Save functional information (eggNOG annotations)
        save_func_data(r, results_store, details_store)

        # Save covariate information (corncob annotations)
        save_stat_data(r, results_store, details_store)

    logger.info("Finished processing")


def copy_table(r, r_key, store, store_key):
    """Read a table from an HDF store and save it to a redis store."""

    logger.info(f"Copying data from {store_key} to {r_key}")

    try:
        r.set(
            r_key,
            pd.read_hdf(
                store,
                store_key
            )
        )
    except KeyError:
        # The store_key does not exist
        logger.info(f"{store_key} does not exist, skipping")
        return


def save_cag_data(r, results_store, details_store):
    """Save information related to CAG composition and abundance."""

    # Save the table with CAG sizes
    copy_table(r, "cag_summary_metrics", results_store, "/annot/cag/all")


def save_expt_data(r, results_store, details_store):
    """Save information related to experimental design."""
    
    # Save the manifest
    r.set("manifest", read_manifest(results_store))

    # Save the number of reads per specimen
    copy_table(r, "specimen_nreads", results_store, "/summary/all")

    # Read the abundance of each CAG across all specimens
    cag_abund_df = pd.read_hdf(
        results_store, 
        "/abund/cag/wide"
    ).set_index(
        "specimen"
    ).drop(
        columns="UNASSIGNED"
    )

    # Save the mean abundance across all specimens
    r.set(
        "mean_abundance_cags",
        cag_abund_df.mean().sort_values(ascending=False)
    )

    # Save the abundance of every CAG across every specimen
    save_sparse_abund(
        r,
        cag_abund_df,
        "cag_abundance_specimen",
        "specimen_abundance_cag"
    )

    # Save the ordination
    r.set(
        "specimen_ordination_pca cag",
        pd.read_hdf(results_store, "/ordination/cag/pca").set_index("specimen")
    )

def save_sparse_abund(r, abund_df, group_key, specimen_key):
    """Save the abundances of every gene group across every specimen in sparse dicts."""

    # For every specimen
    for specimen, c in abund_df.items():

        # Save the abundance of every gene group
        r.set(
            f"{group_key} {specimen}", 
            c.loc[c > 0].to_dict()
        )

    # For every gene group
    for group_id, group_r in abund_df.iterrows():

        # Save the abundance across every specimen
        r.set(
            f"{specimen_key} {group_id}", 
            group_r.loc[group_r > 0].to_dict()
        )


def read_manifest(store):

    logger.info("Reading the manifest")

    # Read the manifest
    manifest = pd.read_hdf(
        store, 
        "/manifest"
    )
    
    # Remove the R1 and R2, if present
    for c in ["R1", "R2"]:
        if c in manifest.columns.values:
            manifest = manifest.drop(columns=[c])

    # Return a deduplicated list    
    return manifest.drop_duplicates(
    ).reset_index(
        drop=True
    )


def save_tax_data(r, results_store, details_store):
    """Save information related to taxonomic annotation of genes."""

    # Read the taxonomic annotation for each gene
    logger.info("Reading in taxonomic assignments")
    try:
        tax_df = pd.read_hdf(results_store, "/annot/gene/tax")
    
    # If there is no taxonomic annotation in this dataset
    except KeyError:
        logger.info("No taxonomic assignments found, skipping")
    
        # then stop any further processing
        return

    # Read the taxonomy
    taxonomy = pd.read_hdf(results_store, "/ref/taxonomy")
    # Format the taxonomy
    taxonomy = taxonomy.apply(
        lambda c: c.apply(int) if c.name == "tax_id" else c
    ).set_index(
        "tax_id"
    )
    # Save the taxonomy
    r.set("taxonomy_df", taxonomy)
    r.set("taxonomy_name", taxonomy["name"].to_dict())
    r.set("taxonomy_rank", taxonomy["rank"].to_dict())

    # Read the abundance of each CAG across all specimens
    taxa_abund_df = pd.read_hdf(
        results_store, 
        "/abund/taxa/wide"
    ).set_index(
        "specimen"
    ).drop(
        columns="UNASSIGNED"
    )

    # Save the mean abundance across all specimens
    r.set(
        "mean_abundance_taxa",
        taxa_abund_df.mean().sort_values(ascending=False)
    )

    # Save the abundance of every taxon across every specimen
    save_sparse_abund(
        r,
        taxa_abund_df,
        "tax_abundance_specimen",
        "specimen_abundance_tax"
    )

    # Save the ordination
    r.set(
        "specimen_ordination_pca taxa",
        pd.read_hdf(results_store, "/ordination/taxa/pca").set_index("specimen")
    )

    # Set the index by gene
    tax_df.set_index("gene", inplace=True)

    # Read the assignment of genes to CAGs
    gene_cag_df = pd.read_hdf(
        results_store, 
        "/annot/gene/cag"
    ).set_index(
        "gene"
    )

    # Get the number of genes assigned to each CAG
    cag_size = gene_cag_df["CAG"].value_counts()

    # Create a Taxonomy object
    logger.info("Creating taxonomy object")
    tax = Taxonomy(results_store)

    # Make a combined DF with the genes that have both
    df = pd.DataFrame(dict(
        CAG=gene_cag_df["CAG"],
        tax_id=tax_df["tax_id"]
    )).dropna(
    ).applymap(
        int
    )
    logger.info(f"Number of genes with both taxonomic and CAG annotations: {df.shape[0]:,}")

    # If no genes have both
    if df.shape[0] == 0:

        # Stop any further processing
        logger.info("No genes have both taxonomic and CAG annotations, skipping")
        return

    # Keep track of the majority assignment at each taxonomic level
    # First key is rank, second key is CAG
    cag_majority_tax = defaultdict(dict)

    # Keep track of the number of genes assigned to each taxon, in a dict
    # keyed by cag ID. This will be used to drive the taxonomic layout of CAGs
    cag_tax_spectra = dict()

    # For each CAG, save a DF with a summary of all taxonomic assignments
    for cag_id, cag_df in df.groupby("CAG"):

        cag_id = int(cag_id)

        # Use the Taxonomy object to format the taxonomic assignment DF
        cag_tax_df = tax.make_cag_tax_df(
            cag_df["tax_id"].value_counts()
        ).sort_values(
            by="total",
            ascending=False
        )

        # Save the table
        r.set(
            f"cag_tax_assignments {cag_id}", 
            cag_tax_df
        )

        # Save the taxonomic spectra
        cag_tax_spectra[
            cag_id
        ] = cag_tax_df.set_index(
            "tax_id"
        )[
            "count"
        ].to_dict()

        # For each taxonomic rank
        for tax_rank, tax_rank_df in cag_tax_df.groupby("rank"):

            # If the top taxon was assigned the majority of genes
            if tax_rank_df["total"].values[0] >= (cag_size[cag_id] / 2):

                # Save this taxon as the majority taxon
                cag_majority_tax[
                    tax_rank
                ][
                    cag_id
                ] = tax_rank_df["name"].values[0]

    # For each rank
    for tax_rank, cag_majority in cag_majority_tax.items():

        # Save the dict of the majority assignment per CAG
        r.set(
            f"cag_majority_taxon {tax_rank}",
            cag_majority
        )

    # For every tax_id
    for tax_id, tax_df in df.groupby("tax_id"):

        # Save a dict with the number of genes from each CAG assigned to this taxon
        r.set(
            f"tax_cag_assignments {tax_id}",
            tax_df["CAG"].value_counts().to_dict()
        )

    # Process the tax spectra
    r.set(
        "cag_taxonomic_layout",
        layout_cags_by_taxonomy(results_store, cag_tax_spectra)
    )
    
def layout_cags_by_taxonomy(results_store, cag_tax_spectra, metric="euclidean", method="complete"):
    """Use the spectrum of taxonomic assignments to perform hierarchical clustering on CAGs."""

    # Get the size of each CAG
    cag_size = pd.read_hdf(results_store, "/annot/cag/all").set_index("cag")["size"]

    # For all CAGs
    for cag_id in cag_size.keys():

        # If there are no genes with taxonomic assignment
        if cag_tax_spectra.get(cag_id) is None:

            # Then add a single value (1) to the 'unassigned' key (-1)
            # This will have the effect of grouping all unassigned CAGs to a single rake
            cag_tax_spectra[cag_id] = {-1: 1}

    # Make a DataFrame with the number of genes assigned to each taxa, per CAG
    logger.info("Computing taxonomic assignments per CAG")
    cag_tax_spectra_df = pd.DataFrame(
        cag_tax_spectra
    ).fillna(
        0
    )

    # Divide by the max to get the proportion
    cag_tax_spectra_df = cag_tax_spectra_df / cag_tax_spectra_df.max()

    # Rotate so that the CAGs are in the index
    cag_tax_spectra_df = cag_tax_spectra_df.T

    # Perform linkage clustering on CAGs (columns)
    logger.info("Performing linkage clustering on CAGs by taxonomic assignment")
    Z = linkage(
        cag_tax_spectra_df,
        method=method,
        metric=metric,
    )

    # Compute the optimal leaf ordering
    Z = optimal_leaf_ordering(
        Z,
        cag_tax_spectra_df,
        metric
    )

    # Use that linkage matrix to build a set of coordinates
    logger.info("Mapping linkage clusters to x-y coordinates")
    pm = PartitionMap(Z, cag_tax_spectra_df.index.values, cag_size)
    coords = pm.get_coords()

    return coords

        
class PartitionMap:
    
    def __init__(self, Z, leaf_names, size_dict):
        """Input is a linkage matrix, and a list of the names of all leaves"""
        
        # Save the dictionary with the number of genes in each subnetwork
        self.size_dict = size_dict

        # Make sure that there is a size value for all listed leaves
        n_missing = np.sum(list(map(lambda n: size_dict.get(n) is None, leaf_names)))
        assert n_missing == 0, f"Missing {n_missing:,}/{len(leaf_names):,} node sizes"

        logger.info(f"Number of leaves: {len(leaf_names):,}")
        
        # Make a DataFrame with the nodes which were joined at each iteration
        Z = pd.DataFrame(
            Z,
            columns = ["childA", "childB", "distance", "size"]
        )

        # Save the number of leaves as a shortcut
        self.nleaves = len(leaf_names)

        # Number internal nodes starting from `nleaves`
        Z = Z.assign(
            node_ix = list(map(
                lambda i: i + self.nleaves,
                np.arange(Z.shape[0])
            ))
        ).apply(
            lambda c: c.apply(int if c.name != 'distance' else float)
        ).set_index(
            "node_ix"
        )

        # Function to get the size of a node / leaf
        # based on the indexing system used by `linkage`
        self.node_size = lambda i: Z.loc[int(i), 'size'] if i >= self.nleaves else size_dict[leaf_names[int(i)]]
        
        # Use the number of genes per subnetwork to update the size of each node
        for node_ix, r in Z.iterrows():

            # Update the table
            try:
                Z.loc[node_ix, 'size'] = self.node_size(r.childA) + self.node_size(r.childB)
            except:
                logger.info(f"Problem updating size: {node_ix} / {r.childA} / {r.childB}")
                logger.info(f"Unexpected error: {sys.exc_info()[0]}")
                raise

        # Reverse the order
        Z = Z[::-1]
        
        # Create the map, assigning the entire range to the final node which was created
        root_node = Z.index.values[0]
        logger.info(f"Root node {root_node} contains {self.node_size(root_node):,} genes")
        self.partition_map = {
            root_node: Partition(root_node, self.node_size(root_node))
        }

        # Save the entire linkage matrix
        self.Z = Z
        
        # Save the leaf names
        self.leaf_names = leaf_names
        
        # Expand each node in turn
        for parent_node_ix, r in Z.iterrows():
            self.add_children(parent_node_ix, r)
            
    def annotate(self, node_ix):
        
        return {
            "size": self.node_size(node_ix),
            "name": node_ix,
            "is_leaf": node_ix < self.nleaves
        }
            
        
    def add_children(self, parent_node_ix, r):
        
        # Make sure that the parent node has a partition assigned
        assert parent_node_ix in self.partition_map
        
        # Check if both children are leaves
        # If they are not leaves, then their size is encoded in the linkage matrix
        childA = self.annotate(r.childA)
        childB = self.annotate(r.childB)
        
        # Split the parent node, weighted by size
        nodeA, nodeB = self.partition_map[
            parent_node_ix
        ].split(
            childA,
            childB,
            r.distance
        )
        
        self.partition_map[childA['name']] = nodeA
        self.partition_map[childB['name']] = nodeB
        
    def get_coords(self):
        """Return the x-y coordinates of all leaf nodes."""
        
        return pd.DataFrame(
            [
                {
                    'name': node.name,
                    'size': node.size,
                    'is_leaf': node.is_leaf,
                    'parent': node.parent,
                    **node.pos()
                }
                for node_name, node in self.partition_map.items()
            ]
        ).astype(
            dict(
                name=int,
                size=int,
                is_leaf=bool,
                parent=int,
                x=float,
                y=float,
            )
        ).set_index(
            "name"
        )



class Partition:

    def __init__(
        self, 
        name,
        size,
        theta_min=0, 
        theta_max=2*np.pi, 
        radius=0, 
        is_leaf=False,
        parent=-1
    ):
        self.name = name
        self.size = size
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.radius = radius
        self.is_leaf = is_leaf
        self.parent = parent
        
    def split(self, childA, childB, radius_delta):
        """Split a node, weighted by child size."""
        # Make sure that the child size adds up to the node size
        assert childA['size'] + childB['size'] == self.size, (childA['size'], childB['size'], self.size)
        assert childA['size'] > 0
        assert childB['size'] > 0
        
        # Get the arc length of the parent
        theta_range = self.theta_max - self.theta_min
        
        # Compute the proportion of that arc which is assigned to each child
        thetaA = theta_range * childA['size'] / self.size
        thetaB = theta_range * childB['size'] / self.size
        
        nodeA = Partition(
            childA['name'],
            childA['size'],
            theta_min = self.theta_min,
            theta_max = self.theta_min + thetaA,
            radius = self.radius + radius_delta,
            is_leaf = childA['is_leaf'],
            parent = self.name
        )

        nodeB = Partition(
            childB['name'],
            childB['size'],
            theta_min = self.theta_max - thetaB,
            theta_max = self.theta_max,
            radius = self.radius + radius_delta,
            is_leaf = childB['is_leaf'],
            parent = self.name
        )

        return nodeA, nodeB
    
    def pos(self):
        """Calculate the position of the node in cartesian coordinates."""
        theta = np.mean([self.theta_min, self.theta_max])
        
        return {
            "name": self.name,
            "x": self.radius * math.cos(theta),
            "y": self.radius * math.sin(theta),
        }

def save_func_data(r, results_store, details_store):
    """Save information related to functional annotation of genes."""

    # Read the functional annotation for each gene
    logger.info("Reading in functional annotations")
    try:
        func_df = pd.read_hdf(results_store, "/annot/gene/eggnog")
    
    # If there is no functional annotation in this dataset
    except KeyError:
        logger.info("No functional annotations found, skipping")
    
        # then stop any further processing
        return

    # Save a dict mapping each `eggnog_desc_ix` integer to an `eggnog_desc` string
    r.set(
        "func_index",
        pd.read_hdf(
            results_store,
            "/ref/eggnog_desc"
        ).set_index(
            "eggnog_desc_ix"
        )[
            "eggnog_desc"
        ].to_dict()
    )

    # Read the abundance of each functional group across all specimens
    func_abund_df = pd.read_hdf(
        results_store, 
        "/abund/func/wide"
    ).set_index(
        "specimen"
    ).drop(
        columns="UNASSIGNED"
    )

    # Save the mean abundance across all specimens
    r.set(
        "mean_abundance_func",
        func_abund_df.mean().sort_values(ascending=False)
    )

    # Save the abundance of every functional group across every specimen
    save_sparse_abund(
        r,
        func_abund_df,
        "specimen_abundance_func",
        "func_abundance_specimen"
    )

    # Save the ordination
    r.set(
        "specimen_ordination_pca func",
        pd.read_hdf(results_store, "/ordination/func/pca").set_index("specimen")
    )

    # Set the index by gene
    func_df.set_index("query_name", inplace=True)

    # Read the assignment of genes to CAGs
    gene_cag_df = pd.read_hdf(
        results_store, 
        "/annot/gene/cag"
    ).set_index(
        "gene"
    )

    # Make a combined DF with the genes that have both
    df = pd.DataFrame(dict(
        CAG=gene_cag_df["CAG"],
        func=func_df["eggnog_desc_ix"]
    )).dropna(
    ).applymap(
        int
    )
    logger.info(f"Number of genes with both functional and CAG annotations: {df.shape[0]:,}")

    # If no genes have both
    if df.shape[0] == 0:

        # Stop any further processing
        logger.info("No genes have both functional and CAG annotations, skipping")
        return

    # For each CAG
    for cag_id, cag_df in df.groupby("CAG"):

        # Save the number of assignments to each function
        r.set(
            f"cag_func_assignments {cag_id}",
            cag_df["func"].value_counts().to_dict()
        )

    # For each function
    for func_ix, func_df in df.groupby("func"):

        # Save the set of CAGs which contain any gene with this assignment
        r.set(
            f"func_cag_set {func_ix}",
            set(func_df["CAG"].tolist())
        )


def save_stat_data(r, results_store, details_store):
    """Save information related to estimation of coefficients of association (corncob)."""

    # Save a list of all parameters
    available_parameters = set()

    # Save stat data for each type of gene grouping
    for group_name in ['cag', 'taxa', 'func']:

        # Read the corncob output for each gene group
        logger.info(f"Reading in corncob results for {group_name}")
        try:
            stat_df = pd.read_hdf(results_store, f"/stats/{group_name}/corncob")

        # If there is no corncob output in this dataset
        except KeyError:
            logger.info(f"No corncob output found for {group_name}, skipping")
        
            # then try the next gene group
            continue

        # Skip the "(Intercept)"
        stat_df = stat_df.query(
            "parameter != '(Intercept)'"
        )

        # For each parameter
        for parameter, parameter_df in stat_df.groupby("parameter"):

            # Add to the list of all available parameters
            available_parameters.add(parameter)

            # Save to redis
            r.set(
                f"{group_name}_association {parameter}",
                parameter_df.drop(
                    columns="parameter"
                ).set_index(
                    group_name
                )
            )

    # Save the list of all parameters
    r.set(
        "available_parameters",
        list(available_parameters)
    )

class Taxonomy:

    def __init__(self, store):
        """Read in the taxonomy table."""

        # Read the taxonomy table
        logger.info("Reading in /ref/taxonomy")

        self.taxonomy_df = pd.read_hdf(
            store,
            "/ref/taxonomy"
        ).apply(
            lambda c: c.fillna(0).apply(float).apply(
                int) if c.name in ["parent", "tax_id"] else c,
        ).set_index(
            "tax_id"
        )

        self.all_taxids = set(self.taxonomy_df.index.values)

    @lru_cache(maxsize=None)
    def path_to_root(self, tax_id, max_steps=100):
        """Parse the taxonomy to yield a list with all of the taxa above this one."""

        visited = []

        for _ in range(max_steps):

            # Skip taxa which are missing
            if tax_id not in self.all_taxids:
                break

            # Add to the path we have visited
            visited.append(tax_id)

            # Get the parent of this taxon
            parent_id = self.taxonomy_df.loc[tax_id, "parent"]

            # If the chain has ended, stop
            if parent_id in visited or parent_id == 0:
                break

            # Otherwise, keep walking up
            tax_id = parent_id

        return visited

    @lru_cache(maxsize=None)
    def anc_at_rank(self, tax_id, rank):
        for anc_tax_id in self.path_to_root(tax_id):
            if self.taxonomy_df.loc[anc_tax_id, "rank"] == rank:
                return anc_tax_id

    def name(self, tax_id):
        if tax_id in self.all_taxids:
            return self.taxonomy_df.loc[tax_id, "name"]

    def parent(self, tax_id):
        if tax_id in self.all_taxids:
            return self.taxonomy_df.loc[tax_id, "parent"]

    def make_cag_tax_df(
        self,
        taxa_vc,
    ):
        """Return a nicely formatted taxonomy table from a list of tax IDs and the number of assignments for each."""

        # We will construct a table with all of the taxa in the tree, containing
        # The ID of that taxon
        # The name of that taxon
        # The number of genes assigned to that taxon (or its children)
        # The rank of that taxon
        # The ID of the parent of that taxon

        # The number of genes found at that taxon or in its decendents
        counts = defaultdict(int)

        # Keep track of the total number of genes with a valid tax ID
        total_genes_assigned = 0

        # Iterate over each terminal leaf
        for tax_id, n_genes in taxa_vc.items():

            # Skip taxa which aren't in the taxonomy
            if tax_id not in self.taxonomy_df.index.values:
                continue

            # Count all genes part of this analysis
            total_genes_assigned += n_genes

            # Walk up the tree from the leaf to the root
            for anc_tax_id in self.path_to_root(tax_id):

                # Add to the sum for every node we visit along the way
                counts[anc_tax_id] += n_genes

        if len(counts) == 0:
            return pd.DataFrame([{
                "tax_id": None
            }])

        # Make a DataFrame
        df = pd.DataFrame({
            "count": counts,
        })

        # Add the name, parent, rank
        df = df.assign(
            tax_id=df.index.values,
            parent=self.taxonomy_df["parent"],
            rank=self.taxonomy_df["rank"],
            name=self.taxonomy_df["name"],
            total=total_genes_assigned,
        ).reset_index(
            drop=True
        )

        return df


# Executed as a script
if __name__ == "__main__":

    # Entrypoint
    buildRedis(
        results=args.results,
        details=args.details,
        host=args.host,
        port=args.port,
    )
