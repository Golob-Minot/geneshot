#!/usr/bin/env python3

import argparse
from collections import defaultdict
from functools import lru_cache
import logging
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as mpatches
import numpy as np
import os
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import cosine
import seaborn as sns
from statsmodels.stats.multitest import multipletests

sns.set_style("whitegrid")


# TAXONOMY OBJECT #
class Taxonomy:
    def __init__(self, hdf_fp):
        self.taxonomy_df = pd.read_hdf(
            hdf_fp,
            "/ref/taxonomy"
        ).apply(
            lambda c: c.apply(str) if c.name == "tax_id" else c
        ).set_index(
            "tax_id"
        )
        self.taxonomy_df["parent"] = self.taxonomy_df["parent"].fillna(
            0).apply(float).apply(int).apply(str)

        # Make a dictionary with the tax IDs for unique names
        name_vc = self.taxonomy_df["name"].value_counts()
        self.name_dict = self.taxonomy_df.loc[
            self.taxonomy_df["name"].isin(name_vc.index.values[name_vc == 1])
        ].reset_index(
        ).set_index(
            "name"
        )["tax_id"]

    def parent(self, tax_id):
        return self.taxonomy_df["parent"].get(str(tax_id))

    def ancestors_iter(self, tax_id, a=[]):

        tax_id = str(tax_id)

        a = a + [tax_id]

        if self.parent(tax_id) is None:
            return a
        elif self.parent(tax_id) == tax_id:
            return a
        else:
            return self.ancestors_iter(
                self.parent(tax_id),
                a=a.copy()
            )

    def name(self, tax_id):
        return self.taxonomy_df["name"].get(str(tax_id))

    def rank(self, tax_id):
        return self.taxonomy_df["rank"].get(str(tax_id))

    def taxid(self, name):
        return self.name_dict.get(name)

    @lru_cache(maxsize=1024)
    def ancestors(self, tax_id):
        return self.ancestors_iter(str(tax_id))

    @lru_cache(maxsize=1024)
    def anc_at_rank(self, tax_id, rank):
        for t in self.ancestors(str(tax_id)):
            if self.rank(t) == rank:
                return t

# Functions for reading taxonomic information
@lru_cache(maxsize=128)
def calc_tax_assignments(cag_id, hdf_fp, tax, max_n_genes=1000):
    # Read in the table with all of the assignments for this CAG
    gene_annots = pd.read_hdf(
        hdf_fp,
        "/annot/gene/all",
        where="CAG == {}".format(cag_id),
        columns=["tax_id"]
    )

    gene_annots = gene_annots.query(
        "tax_id != '0'"
    )

    # If there are more than `max_n_genes` genes, subsample
    if gene_annots.shape[0] > max_n_genes:
        gene_annots = gene_annots.sample(max_n_genes)

    n_genes = gene_annots.shape[0]

    # Count up the number of hits for each tax_id
    tax_id_hits = gene_annots["tax_id"].value_counts()

    # Get the complete set of ancestors which need to be considered
    all_ancestors = list(set([
        a
        for t in tax_id_hits.index.values
        for a in tax.ancestors(t)
        if a != "0"
    ]))

    # Now count up the number of hits which are consistent with each ancestor
    consistent_counts = defaultdict(int)
    # As well as the number of hits which are at-or-below a given ancestor
    at_or_below_counts = defaultdict(int)

    for t, c in tax_id_hits.items():
        for a in all_ancestors:
            if a in tax.ancestors(t):

                consistent_counts[a] += c
                at_or_below_counts[a] += c

            elif t in tax.ancestors(a):
                consistent_counts[a] += c

    df = pd.DataFrame({
        "Specific": at_or_below_counts,
        "Consistent": consistent_counts,
        "Exact": tax_id_hits
    })

    df = df.assign(
        name=list(map(tax.name, df.index.values))
    ).assign(
        rank=list(map(tax.rank, df.index.values))
    ).reset_index(
    ).rename(columns={"index": "tax_id"}).assign(
        Total=n_genes
    )
    df = df.assign(
        prop_consistent=df["Consistent"] / n_genes,
        prop_specific=df["Specific"] / n_genes,
        prop_exact=(df["Exact"] / n_genes).fillna(0),
    )

    return df.sort_values(by="prop_consistent", ascending=False)


# Helper functions to manipulate the taxonomic tree display
def add_taxon(tax_id, terminal_nodes, included_nodes, tax):

    # If we already have this node, skip it
    if tax_id in included_nodes:
        return

    # Add this as a terminal node, and an included node
    terminal_nodes.add(tax_id)
    included_nodes.add(tax_id)

    # Walk through each of the ancestors
    for anc_tax_id in tax.ancestors(tax_id):
        # Only consider ancestors of this taxon
        if anc_tax_id == tax_id:
            continue

        # These higher-level taxa are no longer terminal nodes
        if anc_tax_id in terminal_nodes:
            terminal_nodes.remove(anc_tax_id)
        # All higher level nodes are definitely included
        included_nodes.add(anc_tax_id)


# Function to read from HDF
@lru_cache(maxsize=128)
def read_from_hdf(hdf_fp, hdf_key):
    logging.info("Reading in {} from {}".format(hdf_key, hdf_fp))
    df = pd.read_hdf(hdf_fp, hdf_key)
    if "wald" in df.columns.values:
        df = df.assign(
            abs_wald=df["wald"].abs()
        )
    return df


@lru_cache(maxsize=128)
def calc_functional_assignments(cag_id, hdf_fp):
    # Read in the table with all of the assignments for this CAG
    return pd.read_hdf(
        hdf_fp,
        "/annot/gene/all",
        where="CAG == {}".format(cag_id),
        columns=["eggNOG_desc"]
    ).dropna(
    )["eggNOG_desc"].value_counts()


# Function to select the CAGs to display
@lru_cache(maxsize=1)
def pick_cags_to_plot(
    hdf_fp,
    parameter,
    plot_n_cags=50,
    min_cag_size=10,
    max_cag_size=None
):
    # Now read in the CAG-level corncob results in order to pick the CAGs to plot
    logging.info("Selecting CAGs to plot")
    cag_annot = read_from_hdf(
        hdf_fp,
        "/stats/cag/corncob"
    ).query(
        "parameter == '{}'".format(parameter)
    )
    if "wald" not in cag_annot.columns.values:
        cag_annot = cag_annot.assign(
            wald=cag_annot["estimate"] / cag_annot["std_error"]
        )
        cag_annot = cag_annot.assign(
            abs_wald=cag_annot["wald"].abs()
        )

    # Set the index by CAG ID and add the annotation of CAG size
    cag_annot = cag_annot.set_index(
        "CAG"
    ).assign(
        size=read_from_hdf(
            hdf_fp,
            "/annot/cag/all"
        ).set_index(
            "CAG"
        )[
            "size"
        ]
    ).query(
        "size >= {}".format(min_cag_size)
    )

    if max_cag_size is not None:
        assert isinstance(max_cag_size, int)
        assert max_cag_size > min_cag_size
        cag_annot = cag_annot.query(
            "size <= {}".format(max_cag_size)
        )
    
    cag_annot = cag_annot.sort_values(
        by="abs_wald",
        ascending=False
    ).head(
        plot_n_cags
    )
    cags_to_plot = cag_annot.index.values

    return cags_to_plot, cag_annot


# Function to select the eggNOG functions to display
def get_betta_functions(
    hdf_fp,
    parameter,
    plot_n_functions=35,
    plot_n_cags=50,
    min_cag_size=10,
    exclude_prefixes=[],
    exclude_terms=[],
    show_only_functions=None,
    cags_to_plot=None,
    max_function_string_length=70
):
    # Read in the betta results
    # Filter by parameter
    # Sort by absolute wald
    # Filter out everything except the eggNOG results
    betta = read_betta(
        hdf_fp, 
        parameter
    ).query(
        "annotation == 'eggNOG_desc'"
    ).set_index(
        "label"
    ).sort_values(
        by="abs_wald",
        ascending=False
    )
    assert betta.shape[0] > 0

    betta_wald = betta["wald"]

    # Only pick the CAGs to plot if they haven't been otherwise provided
    if cags_to_plot is None:
        # If the CAGs have not been provided externally, cluster the plot by linkage
        col_cluster = True
        # Select the top `top_n_cags` CAGs with the highest abs(wald), but filtering on CAG size
        cags_to_plot, _ = pick_cags_to_plot(
            hdf_fp,
            parameter,
            plot_n_cags=plot_n_cags,
            min_cag_size=min_cag_size,
        )
    else:
        # If the CAGs have been provided, don't cluster the columns
        col_cluster = False

    # Get the functional annotations for those cags
    logging.info("Reading functional assignments for CAGs")
    cag_func_df = pd.DataFrame({
        cag_id: calc_functional_assignments(cag_id, hdf_fp)
        for cag_id in cags_to_plot
    }).fillna(0).reindex(
        columns=cags_to_plot
    )

    # If show_only_functions was specified, use that as the list of functions to plot
    if show_only_functions is not None:
        funcs_to_plot = [
            n for n in show_only_functions.split(",")
            if len(n) > 0
        ]
        assert len(funcs_to_plot) > 0, "Please provide functions with --show-only-functions (none found)"

        for func_name in funcs_to_plot:
            assert func_name in cag_func_df.index.values, "{} not found in results".format(func_name)

    else:
        # Otherwise, pick the top functions to plot as the top functions by abs_wald as a subset of those in these CAGs
        funcs_to_plot = []
        for d in betta_wald.reindex(
            index=cag_func_df.index.values
        ).dropna().abs().sort_values(ascending=False).index.values:
            # Skip excluded terms
            if len(exclude_terms) > 0 and d in exclude_terms:
                continue
            # Skip excluded prefixes
            if any([d.startswith(n) for n in exclude_prefixes]):
                continue
            # Stop when the list is full
            if len(funcs_to_plot) >= plot_n_functions:
                break
            # Add to the list
            if d not in funcs_to_plot:
                funcs_to_plot.append(d)

    print("\n\n".join(funcs_to_plot))

    # Subset to these functions
    cag_func_df = cag_func_df.reindex(
        index=funcs_to_plot,
    ).fillna(0)

    # Transform to boolean
    cag_func_df = cag_func_df > 0

    # Filter the functions down by string length
    cag_func_df = cag_func_df.rename(
        index=lambda n: n[:max_function_string_length] if len(
            n) > max_function_string_length else n
    )

    # Reorder the functions by co-occurance
    cag_func_df = cag_func_df.reindex(
        index=cag_func_df.index.values[
            leaves_list(linkage(
                cag_func_df,
                method="ward",
                metric="euclidean",
                optimal_ordering=True,
            ))
        ]
    )

    # Optionally reorder the columns (CAGs)
    if col_cluster:
        cag_func_df = cag_func_df.reindex(
            columns=cag_func_df.columns.values[
                leaves_list(linkage(
                    cag_func_df.T,
                    method="ward",
                    metric="euclidean",
                    optimal_ordering=True,
                ))
            ]
        )

    return cag_func_df, betta.reindex(index=cag_func_df.index, columns=["estimate", "std_error", "wald"])


def read_betta(hdf_fp, parameter):
    df = read_from_hdf(
        hdf_fp,
        "/stats/enrichment/betta"
    ).query(
        "parameter == '{}'".format(parameter)
    )
    if "abs_wald" not in df.columns.values:
        df = df.assign(
            wald = df["estimate"] / df["std_error"]
        )
        df = df.assign(
            abs_wald = df["wald"].abs()
        )

    return df.sort_values(
        by="abs_wald",
        ascending=False
    )


def plot_betta_taxa_and_functions(
    hdf_fp,
    parameter,
    tax,
    betta,
    cags_to_plot,
    cag_annot,
    pdf=None,
    taxonomic_ranks=["family", "genus", "species"],
    plot_n_taxa=20,
    plot_n_cags=50,
    plot_n_functions=25,
    min_cag_size=20,
    rank_order=["phylum", "class", "order", "family", "genus", "species"],
    figsize=[22, 18],
    exclude_terms=[
        'ATPases associated with a variety of cellular activities',
        'domain protein',
        'Transcriptional regulator',
        'domain, Protein',
        'Required for chromosome condensation and partitioning',
    ],
    exclude_prefixes=[
        "Psort location"
    ],
    include_functions=True,
    show_only_functions=None,
):

    # Filter out the non-taxonomic results
    assert betta["annotation"].isin(taxonomic_ranks).any(), "No annotations found with {}".format(" / ".join(taxonomic_ranks))
    betta = betta.loc[
        betta["annotation"].isin(taxonomic_ranks)
    ]

    # Get the taxonomic annotations for those cags
    logging.info("Reading taxonomic assignments for CAGs")
    cag_taxa_df = pd.concat([
        calc_tax_assignments(cag_id, hdf_fp, tax).assign(
            CAG=cag_id
        )
        for cag_id in cags_to_plot
    ])

    # Use that set of data to assign a tax ID to the betta results, where possible
    logging.info("Tabulating common tax ID mappings")
    taxa_name_dict = cag_taxa_df.reindex(
        columns=["tax_id", "name"]
    ).drop_duplicates(
    ).set_index("name")["tax_id"]

    # Check that dictionary and, as a backup, query the taxonomy for a unique match to the name
    logging.info("Assigning tax IDs to betta results")
    betta = betta.assign(
        tax_id=betta["label"].apply(
            lambda label: taxa_name_dict.get(label, tax.taxid(label))
        )
    )
    tax_dict = betta.set_index("label")["tax_id"]

    # Pick the nodes to display by walking down the list of results
    # and adding as we go, keeping track of how many different terminal nodes are in the plot
    terminal_nodes = set([])
    included_nodes = set([])
    for _, betta_r in betta.iterrows():
        assert betta_r["tax_id"] is not None, "Need to manually add tax ID {}".format(
            betta_r["tax_id"])

        add_taxon(betta_r["tax_id"], terminal_nodes, included_nodes, tax)

        # Stop when we have found the maximum number of terminal nodes
        if len(terminal_nodes) >= plot_n_taxa:
            break

    # Now we need to go through and make sure to add any taxon which is >= 25% of the selected CAGs
    for tax_id in cag_taxa_df.query(
        "prop_specific >= 0.25"
    )["tax_id"].drop_duplicates().values:
        add_taxon(tax_id, terminal_nodes, included_nodes, tax)

    # Make a DataFrame with all of the included taxa for the betta results
    betta_taxa_df = pd.DataFrame({
        tax.name(terminal_tax_id): {
            rank: tax.anc_at_rank(terminal_tax_id, rank)
            for rank in rank_order
        }
        for terminal_tax_id in list(terminal_nodes)
    }).T

    # Sort alphabetically by scientific name
    betta_taxa_df = betta_taxa_df.reindex(
        index=betta_taxa_df.applymap(
            tax.name
        ).sort_values(
            by=rank_order
        ).index
    )

    # To make plotting more understandible, we can fill the most specific taxonomic rank
    # with the name of whatever organism has been measured immediately above it
    betta_taxa_df = betta_taxa_df.apply(
        lambda c: c.fillna(
            betta_taxa_df.apply(
                lambda r: r.dropna(
                ).values[-1] if r.dropna().shape[0] > 0 else "",
                axis=1
            ).apply(
                lambda s: " {} ".format(s)
            )
        ) if c.name == rank_order[-1] else c
    )

    # Drop any rows with no taxa assigned
    betta_taxa_df = betta_taxa_df.reindex(
        index=[
            i
            for i in betta_taxa_df.index.values
            if pd.isnull(i) is False
        ]
    )

    # For each taxon, figure out the x- and y-coordinates from the taxonomy table
    plot_dat = []
    x_pos = 0
    for rank_name, rank_taxa in betta_taxa_df.items():
        x_pos -= 2

        for org_name in rank_taxa.dropna().unique():
            # The index position is the start of the vertical range
            y_start = rank_taxa.tolist().index(org_name)
            # The number of cells with this organism in the height of the vertical range
            y_height = (rank_taxa == org_name).sum()

            # The y position is in the middle of the vertical range
            y_pos = y_start + (y_height / 2)

            # Get the first non-null ancestor for this node
            if rank_order.index(rank_name) == 0:
                ancestor = None
            else:
                ancestor = betta_taxa_df.loc[
                    betta_taxa_df[rank_name] == org_name, :rank_name
                ].iloc[
                    0
                ].dropna()
                if ancestor.shape[0] > 1:
                    ancestor = ancestor.values[-2]
                else:
                    ancestor = None

            plot_dat.append({
                "name": tax.name(org_name.strip(" ")),
                "tax_id": org_name,
                "x": x_pos,
                "y": y_pos,
                "ancestor": ancestor if ancestor != org_name else None,
                "rank": rank_name
            })

    plot_dat = pd.DataFrame(plot_dat)

    # Add the wald statistic
    plot_dat = plot_dat.assign(
        Wald=plot_dat["name"].fillna(
            ""
        ).apply(
            lambda s: s.strip(" ")
        ).apply(
            betta.set_index("label")["wald"].get
        ).fillna(0)
    ).set_index(
        "tax_id"
    )

    # SETTING UP THE PLOT AREA

    # SET UP DIFFERENTLY DEPENDING ON WHETHER THE FUNCTIONS WILL BE INCLUDED
    if include_functions:
        fig, axarr = plt.subplots(
            4,
            3,
            figsize=figsize,
            sharey="row",
            gridspec_kw={
                'height_ratios': [6, 6, 1, 1],
                'width_ratios': [3, 1, 3],
            }
        )

        heatmap_ax = axarr[0, 0]
        taxa_wald_ax = axarr[0, 1]
        tree_ax = axarr[0, 2]
        func_ax = axarr[1, 0]
        func_wald_ax = axarr[1, 1]
        cag_wald_ax = axarr[2, 0]
        cag_size_ax = axarr[3, 0]

        empty_subplots = [
            axarr[1, 2],
            axarr[2, 1],
            axarr[2, 2],
            axarr[3, 1],
            axarr[3, 2],
        ]

    else:
        fig, axarr = plt.subplots(
            3,
            3,
            figsize=figsize,
            sharey="row",
            gridspec_kw={
                'height_ratios': [6, 1, 1],
                'width_ratios': [3, 1, 3],
            }
        )
        heatmap_ax = axarr[0, 0]
        taxa_wald_ax = axarr[0, 1]
        tree_ax = axarr[0, 2]
        cag_wald_ax = axarr[1, 0]
        cag_size_ax = axarr[2, 0]

        empty_subplots = [
            axarr[1, 1],
            axarr[1, 2],
            axarr[2, 1],
            axarr[2, 2],
        ]

    # Clear the unused subplots
    for ax in empty_subplots:
        ax.axis("off")
        ax.grid(False)

    # Add the lines and labels to the taxa plot
    for org_name, r in plot_dat.iterrows():
        # Label each organism
        tree_ax.text(
            r["x"],
            r["y"] + 0.5,
            tax.name(org_name.strip(" ")),
            horizontalalignment="center",
            verticalalignment="center",
            rotation=5,
            zorder=2,
        )

        if r["ancestor"] is not None:

            # Add a line to the ancestor
            if r["x"] - plot_dat.loc[r["ancestor"], "x"] == 1.:
                tree_ax.plot(
                    [r["x"], plot_dat.loc[r["ancestor"], "x"]],
                    [r["y"], plot_dat.loc[r["ancestor"], "y"]],
                    linestyle="--",
                    alpha=0.5,
                    color="grey",
                    zorder=1,
                )
            else:
                tree_ax.plot(
                    [r["x"], plot_dat.loc[r["ancestor"], "x"] -
                        1, plot_dat.loc[r["ancestor"], "x"]],
                    [r["y"], r["y"], plot_dat.loc[r["ancestor"], "y"]],
                    linestyle="--",
                    alpha=0.5,
                    color="grey",
                    zorder=1,
                )

    # Draw the points for the tree over the lines
    sns.scatterplot(
        data=plot_dat,
        x="x",
        y="y",
        hue="Wald",
        hue_norm=(-plot_dat["Wald"].abs().max(), plot_dat["Wald"].abs().max()),
        palette="RdBu",
        linewidth=1,
        edgecolor="black",
        zorder=2,
        ax=tree_ax
    )

    tree_ax.axis("off")
    tree_ax.grid(False)
    tree_ax.legend(bbox_to_anchor=[1.1, 0.9])

    # Now make a plot with the proportion of genes assigned to the top CAGs
    cag_prop_df = pd.DataFrame({
        cag_id: {
            row_ix: d.loc[
                d["name"].isin(row_taxa.apply(
                    tax.name).dropna().drop_duplicates().values),
                "prop_exact"
            ].sum()
            for row_ix, row_taxa in betta_taxa_df.iterrows()
        }
        for cag_id, d in cag_taxa_df.groupby("CAG")
    })

    # Reorder the DataFrame for plotting to match the taxonomic tree
    cag_prop_df = cag_prop_df.reindex(
        index=betta_taxa_df.index.values,
        columns=cag_prop_df.columns.values[
            leaves_list(linkage(
                cag_prop_df.T,
                method="ward",
                metric="euclidean",
                optimal_ordering=True,
            ))
        ]
    )

    # PLOT THE CAG TAXONOMIC HEATMAP
    sns.heatmap(
        data=cag_prop_df.rename(
            columns=lambda cag_id: "CAG {}".format(cag_id)
        ),
        cmap="Blues",
        ax=heatmap_ax,
        cbar=False,
        xticklabels=1,
        yticklabels=1
    )

    heatmap_ax.xaxis.set_ticks_position('top')
    heatmap_ax.tick_params(
        axis="x",
        rotation=90,
    )

    # PLOT THE WALD BY TERMINAL NODE IN TAXA GRAPH
    taxa_wald = betta.set_index("label").reindex(
        index=cag_prop_df.index.values
    )["wald"]

    taxa_wald.plot(
        kind="barh",
        ax=taxa_wald_ax,
        align="edge",
        width=1,
    )
    taxa_wald_ax.xaxis.set_ticks_position('top')
    taxa_wald_ax.set_title("Wald Statistic")
    taxa_wald_ax.tick_params(
        axis="x",
        labelbottom=False,
        bottom=False,
        labeltop=True,
        top=True
    )

    if include_functions:
        # MAKE THE CAG FUNCTION HEATMAP
        betta_func_df, betta_func_details = get_betta_functions(
            hdf_fp,
            parameter,
            exclude_terms=exclude_terms,
            exclude_prefixes=exclude_prefixes,
            show_only_functions=show_only_functions,
            cags_to_plot=cag_prop_df.columns.values,
            plot_n_functions=plot_n_functions
        )

        # PLOT THE CAG FUNCTION WALD STATISTICS
        betta_func_details["wald"].plot(
            kind="barh",
            width=1,
            align="edge",
            ax=func_wald_ax,
        )
        func_wald_ax.set_xlabel("Wald Statistic")

        # Set the xlim to be the same for both of the Wald barplots for taxa and functions
        wald_min = min(betta_func_details["wald"].min(), taxa_wald.min(), 0)
        wald_max = max(betta_func_details["wald"].max(), taxa_wald.max(), 0)
        wald_span = wald_max - wald_min
        func_wald_ax.set_xlim(
            wald_min - (wald_span * 0.05),
            wald_max + (wald_span * 0.05)
        )
        taxa_wald_ax.set_xlim(
            wald_min - (wald_span * 0.05),
            wald_max + (wald_span * 0.05)
        )

        sns.heatmap(
            data=betta_func_df,
            cmap="Blues",
            ax=func_ax,
            cbar=False,
            xticklabels=[],
            yticklabels=1,
        )

        # Rotate the yticklabels
        func_ax.set_yticklabels(func_ax.get_yticklabels(), rotation=0)

    # Make a plot with the CAG size in the lowest subplot
    cag_annot.reindex(
        cag_prop_df.columns.values
    )["size"].rename(
        index=lambda cag_id: "CAG {}".format(cag_id)
    ).apply(
        np.log10
    ).plot(
        kind="bar",
        ax=cag_size_ax,
        width=1,
    )

    # Make a plot with the CAG wald statistic in the second lowest subplot
    cag_annot.reindex(
        cag_prop_df.columns.values
    )["wald"].plot(
        kind="bar",
        ax=cag_wald_ax,
        width=1,
    )

    # Customize the axis labels
    cag_size_ax.set_ylabel(
        "CAG Size     \n(# of genes, log10)",
    )
    cag_wald_ax.set_xticks([])
    cag_wald_ax.set_xlabel("")
    cag_wald_ax.set_ylabel(
        "Wald Statistic",
    )

    # # Set the limits of the horizontal CAG axis
    # cag_wald_ax.set_xlim(-0.5, cag_annot.shape[0] - 0.5)
    # cag_size_ax.set_xlim(-0.5, cag_annot.shape[0] - 0.5)

    plt.tight_layout(w_pad=0.05, h_pad=0.1)

    # Adjust the vertical position of the row labels for the taxonomic heatmap
    heatmap_ax.set_yticks([
        v + 0.5
        for v in range(cag_prop_df.shape[0])
    ])

    if pdf is not None:
        pdf.savefig(bbox_inches="tight")

    plt.close()


def plot_top_annotations(betta, pdf=None, exclude_terms=[], exclude_prefixes=[]):
    for annotation, annotation_df in betta.groupby("annotation"):
        plot_df = annotation_df.copy(
        ).set_index(
            "label"
        )

        # Remove any blacklisted annotations
        to_remove = []
        for d in plot_df.index.values:
            # Skip excluded terms
            if len(exclude_terms) > 0 and d in exclude_terms:
                to_remove.append(d)
                continue
            # Skip excluded prefixes
            if any([d.startswith(n) for n in exclude_prefixes]):
                to_remove.append(d)
                continue
        if len(to_remove) > 0:
            plot_df = plot_df.drop(index=to_remove)

        fig, ax = plt.subplots()
        ax.axvline(x=0, linestyle="--", alpha=0.5)

        y = 0
        label_list = []
        for label, r in plot_df.sort_values(
            by="abs_wald",
            ascending=False
        ).head(
            20
        ).iterrows(
        ):
            ax.plot(
                [r["estimate"] - (r["std_error"] / 2), r["estimate"] + (r["std_error"] / 2)],
                [y, y],
                color="black"
            )
            ax.scatter([r["estimate"]], [y], color="black")
            y += 1
            label_list.append(label)

        ax.set_yticks(list(range(len(label_list))))
        ax.set_yticklabels(label_list)
        ax.set_title(annotation)
        plt.tight_layout()

        if pdf is not None:
            pdf.savefig(bbox_inches="tight")

        plt.close()

        
if __name__ == "__main__":

    log_formatter = logging.Formatter(
        "%(asctime)s %(levelname)-8s [Geneshot Plot Betta] %(message)s"
    )
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    # Write logs to STDOUT
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    parser = argparse.ArgumentParser(
        description="""
        Plot estimated coefficients aggregated by taxa and function using betta.

        Note: compatible with output from geneshot v0.6.0 and higher.

        Example Usage:

        plot_betta_associations.py \
            --hdf <HDF_FP> \
            --parameter <PARAMETER> \
            --out <PDF_FP>

        Required:
            --hdf: Path to results HDF5 file generated by geneshot
            --parameter: Name of the parameter from the formula to be displayed
            --out: Path to output PDF

        Run with --help for a complete list of optional parameters.
        """
    )

    parser.add_argument(
        "--hdf",
        type=str,
        required=True,
        help="Path to results HDF5 file generated by geneshot"
    )

    parser.add_argument(
        "--parameter",
        type=str,
        required=True,
        help="Name of the parameter from the formula to be displayed"
    )

    parser.add_argument(
        "--out",
        type=str,
        required=True,
        help="Path to output PDF"
    )

    parser.add_argument(
        "--cag-annotations",
        type=str,
        default=None,
        help="If specified, write out a table of annotations for all displayed CAGs in CSV format"
    )

    parser.add_argument(
        "--n-cags",
        type=int,
        default=50,
        help="Number of CAGs to include"
    )

    parser.add_argument(
        "--min-cag-size",
        type=int,
        default=10,
        help="Only display CAGs containing at least this many genes"
    )

    parser.add_argument(
        "--max-cag-size",
        type=int,
        default=None,
        help="If specified, exclude CAGs made up of more than this number of genes"
    )

    parser.add_argument(
        "--n-taxa",
        type=int,
        default=20,
        help="Number of taxa to include"
    )

    parser.add_argument(
        "--n-functions",
        type=int,
        default=50,
        help="Number of eggNOG functions to include"
    )

    parser.add_argument(
        "--exclude-prefixes",
        type=str,
        default="Psort location;Required for replicative DNA synthesis;In addition to polymerase activity;DNA-dependent RNA polymerase;Couples transcription and DNA repair;Is required not only for elongation of protein synthesis but also;Protein of unknown function",
        help="Exclude any eggNOG functions starting with these strings (semi-colon delimited)"
    )

    parser.add_argument(
        "--exclude-functions",
        type=str,
        default="ATPases associated with a variety of cellular activities;domain protein;Transcriptional regulator;domain, Protein;Required for chromosome condensation and partitioning;response regulator receiver;Transcriptional regulatory protein, C terminal;Belongs to ABC transporter superfamily;response regulator;Part of a membrane complex involved in electron transport;non supervised orthologous group",
        help="Exclude any eggNOG functions matching these strings (semi-colon delimited)"
    )

    parser.add_argument(
        "--show-only-functions",
        type=str,
        default=None,
        help="If specified, only plot the functions which are explicitly defined in this list (comma-delimited, surrounded by quotes to capture whitespace)."
    )

    parser.add_argument(
        "--figsize",
        type=str,
        default="22,18",
        help="Width and height of plot"
    )

    # Parse the arguments
    args = parser.parse_args()

    logging.info("HDF5 Input File: {}".format(args.hdf))
    logging.info("Parameter: {}".format(args.parameter))
    logging.info("Output PDF: {}".format(args.out))
    
    # Parse the figsize
    figsize = [int(n) for n in args.figsize.split(',')]
    assert len(figsize) == 2

    # Set up the taxonomy
    logging.info("Reading in the taxonomy")
    tax = Taxonomy(args.hdf)

    # Select the top `top_n_cags` CAGs with the highest abs(wald), but filtering on CAG size
    cags_to_plot, cag_annot = pick_cags_to_plot(
        args.hdf,
        args.parameter,
        plot_n_cags=args.n_cags,
        min_cag_size=args.min_cag_size,
        max_cag_size=args.max_cag_size,
    )

    # If specified, write out the annotations for these CAGs
    if args.cag_annotations is not None:
        logging.info(
            "Preparing to write CAG annotations to {} in CSV format".format(
                args.cag_annotations
            )
        )
        pd.concat([
            pd.read_hdf(
                args.hdf,
                "/annot/gene/all",
                where="CAG == {}".format(cag_id),
            )
            for cag_id in cags_to_plot
        ]).to_csv(
            args.cag_annotations,
            index=None
        )
        logging.info("Done writing to {}".format(args.cag_annotations))

    # Read in the betta results
    # Filter by parameter
    # Sort by absolute wald
    betta = read_betta(args.hdf, args.parameter)

    # Make the plot
    with PdfPages(args.out) as pdf:

        plot_top_annotations(
            betta, 
            pdf=pdf,
            exclude_terms=args.exclude_functions.split(";"),
            exclude_prefixes=args.exclude_prefixes.split(";"),
        )

        for include_functions in [False, True]:
            plot_betta_taxa_and_functions(
                args.hdf, 
                args.parameter, 
                tax,
                betta,
                cags_to_plot,
                cag_annot,
                pdf=pdf,
                exclude_terms=args.exclude_functions.split(";"),
                exclude_prefixes=args.exclude_prefixes.split(";"),
                show_only_functions=args.show_only_functions,
                plot_n_taxa=args.n_taxa,
                plot_n_cags=args.n_cags,
                plot_n_functions=args.n_functions,
                taxonomic_ranks=["family", "genus", "species"],
                min_cag_size=args.min_cag_size,
                rank_order = ["phylum", "class", "order", "family", "genus", "species"],
                figsize=[
                    figsize[0],
                    figsize[1] if include_functions else figsize[1] * 0.65
                ],
                include_functions=include_functions
            )
