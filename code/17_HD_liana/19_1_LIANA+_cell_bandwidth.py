#!/usr/bin/env python3

# Try different bandwidths to see the results

import argparse
import json
import os

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotnine as p9
import scanpy as sc
import seaborn as sns
import squidpy as sq

import liana as li
import muon as mu
import decoupler as dc

from pyhere import here


###############################################################################
# Command-line arguments
###############################################################################

parser = argparse.ArgumentParser()

parser.add_argument(
    "--bandwidth",
    type=float,
    required=True
)

args = parser.parse_args()

bandwidth = args.bandwidth

print("\n========================================")
print(f"bandwidth = {bandwidth}")
print("========================================\n")


###############################################################################
# Paths
###############################################################################

in_path = here(
    "processed-data",
    "10_HD_bin_level",
    "no_secondary",
    "liana"
)

out_path = here(
    "processed-data",
    "10_HD_bin_level",
    "no_secondary",
    "liana2"
)

os.makedirs(
    out_path,
    exist_ok=True
)

figure_path = here(
    "plots",
    "10_HD_bin_level",
    "no_secondary",
    "liana2"
)

os.makedirs(
    figure_path,
    exist_ok=True
)


###############################################################################
# Select sample from SLURM array
###############################################################################

sample_info_path = here(
    "raw-data",
    "sample_info",
    "hd_basic_info_split.csv"
)

sample_info = pd.read_csv(
    sample_info_path
)

slurm_task_id = int(
    os.getenv("SLURM_ARRAY_TASK_ID")
)

sample_id = sample_info["sample_id"].iloc[
    slurm_task_id - 1
]

tissue_id = sample_info["tissue_id"].iloc[
    slurm_task_id - 1
]

print(f"SLURM_ARRAY_TASK_ID = {slurm_task_id}")
print(f"sample_id = {sample_id}")
print(f"tissue_id = {tissue_id}")


###############################################################################
# Scale-factor JSON path
###############################################################################

scale_json = here(
    "processed-data",
    "01_spaceranger",
    "five_samples_10_2025",
    sample_id,
    "outs",
    "binned_outputs",
    "square_008um",
    "spatial",
    "scalefactors_json.json"
)


###############################################################################
# Read AnnData
###############################################################################

file_path = os.path.join(
    in_path,
    "adata",
    "adata_cellular_withcelltype_newcluster.h5ad"
)

adata = sc.read(
    file_path
)

adata.var.index.name = None


###############################################################################
# Run separately for each sample
###############################################################################

adata.obs = adata.obs.rename(
    columns={
        "tissue_section": "tissue_id"
    }
)

adata = adata[
    adata.obs["tissue_id"].isin([tissue_id])
].copy()

print(
    f"\nNumber of observations after tissue subset: "
    f"{adata.n_obs}"
)


###############################################################################
# Gene / cell filtering
###############################################################################

# Remove deprecated genes
adata = adata[
    :,
    ~adata.var_names.str.startswith("DEPRECATED_")
].copy()

adata.var_names = adata.var_names.astype(str)

adata.var_names_make_unique()


# Remove bins/cells labelled "Drop"
adata = adata[
    adata.obs["fine_cell_type"] != "Drop"
].copy()

print(
    f"Number of observations after removing Drop: "
    f"{adata.n_obs}"
)


###############################################################################
# Plot fine cell types in spatial coordinates
###############################################################################

ax = sc.pl.embedding(
    adata,
    basis="spatial",
    color=["fine_cell_type"],
    wspace=0.4,
    s=5,
    show=False
)

ax.invert_yaxis()

ax.set_aspect(
    "equal",
    adjustable="box"
)

plt.savefig(
    os.path.join(
        figure_path,
        f"{tissue_id}_fine_cell_type_spatial.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close()


###############################################################################
# Basic preparation and QC
###############################################################################

# Filter cells/bins
sc.pp.filter_cells(
    adata,
    min_genes=10
)

# Filter genes
sc.pp.filter_genes(
    adata,
    min_cells=3
)

print(
    f"\nAfter QC:"
    f"\n  n_obs = {adata.n_obs}"
    f"\n  n_vars = {adata.n_vars}"
)


###############################################################################
# Store counts
###############################################################################

adata.layers["counts"] = adata.X.copy()


###############################################################################
# Normalize expression
###############################################################################

sc.pp.normalize_total(
    adata,
    target_sum=1e4
)

sc.pp.log1p(
    adata
)


###############################################################################
# Check spatial coordinates
###############################################################################

print("\nSpatial coordinate information:")

print(
    f"  coordinate matrix shape = "
    f"{adata.obsm['spatial'].shape}"
)

print(
    f"  x range = "
    f"{adata.obsm['spatial'][:, 0].min():.2f} "
    f"to "
    f"{adata.obsm['spatial'][:, 0].max():.2f}"
)

print(
    f"  y range = "
    f"{adata.obsm['spatial'][:, 1].min():.2f} "
    f"to "
    f"{adata.obsm['spatial'][:, 1].max():.2f}"
)


###############################################################################
# Query bandwidth
###############################################################################

plot, df = li.ut.query_bandwidth(
    coordinates=adata.obsm["spatial"],
    start=0,
    end=adata.n_obs - 1,
    interval_n=40
)


###############################################################################
# Plot bandwidth vs number of neighbours
###############################################################################

ymin = int(
    df.neighbours.min()
)

ymax = int(
    df.neighbours.max()
)

step = max(
    1,
    round(
        (ymax - ymin) / 6
    )
)

p = (
    plot
    + p9.geom_vline(
        xintercept=bandwidth,
        linetype="dashed",
        color="red"
    )
    + p9.annotate(
        "text",
        x=bandwidth,
        y=ymax,
        label=f"chosen = {bandwidth:.0f}",
        angle=90,
        va="bottom",
        ha="right",
        color="red"
    )
    + p9.scale_y_continuous(
        breaks=list(
            range(
                ymin,
                ymax + 1,
                step
            )
        )
    )
    + p9.coord_cartesian(
        xlim=(0, 10000)
    )
)

p.save(
    filename=os.path.join(
        figure_path,
        f"{tissue_id}_spatial_connectivity_bandwidth_{bandwidth}.pdf"
    ),
    dpi=300,
    width=6,
    height=4,
    units="in"
)


###############################################################################
# Calculate spatial connectivity matrix Wij
###############################################################################

li.ut.spatial_neighbors(
    adata=adata,
    bandwidth=bandwidth,
    spatial_key="spatial",
    max_neighbours=adata.n_obs - 1
)


###############################################################################
# Select one spatial point for connectivity visualization
###############################################################################

idx = min(
    5500,
    adata.n_obs - 1
)


###############################################################################
# Get exact coordinates of selected point
###############################################################################

x_coord = adata.obsm["spatial"][idx, 0]
y_coord = adata.obsm["spatial"][idx, 1]

obs_name = adata.obs_names[idx]


###############################################################################
# Get metadata for selected point
###############################################################################

fine_cell_type = (
    adata.obs.iloc[idx]["fine_cell_type"]
    if "fine_cell_type" in adata.obs.columns
    else np.nan
)

region = (
    adata.obs.iloc[idx]["region"]
    if "region" in adata.obs.columns
    else np.nan
)


###############################################################################
# Print selected point information
###############################################################################

print("\n========================================")
print("Selected connectivity center")
print("========================================")

print(f"idx              = {idx}")
print(f"obs_name         = {obs_name}")
print(f"x coordinate     = {x_coord:.4f}")
print(f"y coordinate     = {y_coord:.4f}")
print(f"fine_cell_type   = {fine_cell_type}")
print(f"region           = {region}")
print(f"tissue_id        = {tissue_id}")
print(f"bandwidth        = {bandwidth}")

print("========================================\n")


###############################################################################
# Save selected point information
###############################################################################

selected_point_df = pd.DataFrame(
    {
        "tissue_id": [tissue_id],
        "sample_id": [sample_id],
        "bandwidth": [bandwidth],
        "idx": [idx],
        "obs_name": [obs_name],
        "x_coordinate": [x_coord],
        "y_coordinate": [y_coord],
        "fine_cell_type": [fine_cell_type],
        "region": [region]
    }
)

selected_point_file = os.path.join(
    out_path,
    f"{tissue_id}_connectivity_center_{bandwidth}.csv"
)

selected_point_df.to_csv(
    selected_point_file,
    index=False
)

print(
    f"Saved selected point information:\n"
    f"{selected_point_file}"
)


###############################################################################
# Visualize proximity/connectivity graph
###############################################################################

p = (
    li.pl.connectivity(
        adata,
        idx=idx,
        size=0.01,
        spatial_key="spatial",
        return_fig=True
    )
    + p9.coord_fixed()
)

p.save(
    filename=os.path.join(
        figure_path,
        f"{tissue_id}_connectivity_idx{idx}_{bandwidth}.pdf"
    ),
    dpi=300,
    width=6,
    height=10,
    units="in"
)


###############################################################################
# Spatially variable genes
###############################################################################

# Current analysis DOES NOT restrict LIANA analysis to SVGs.
#
# Example Moran's I approach:
#
# sq.gr.spatial_autocorr(
#     adata,
#     mode="moran",
#     use_raw=False,
#     show_progress_bar=True
# )
#
# svgs = adata.uns["moranI"].index[
#     (
#         adata.uns["moranI"]["pval_norm_fdr_bh"] < 0.05
#     )
#     &
#     (
#         adata.uns["moranI"]["I"] > 0.01
#     )
# ]
#
# print(len(svgs))
#
# adata = adata[:, svgs]


# Alternative: previously generated SVG list
#
# list_path = here(
#     "processed-data",
#     "10_HD_bin_level",
#     "no_secondary",
#     "nnSVG_out",
#     "merged_SVGs.txt"
# )
#
# with open(list_path, "r") as f:
#     svgs = f.read().splitlines()
#
# svg_mask = adata.var["gene_id"].isin(svgs)
#
# adata = adata[:, svg_mask]


###############################################################################
# LIANA inflow analysis
###############################################################################

resource = li.rs.select_resource(
    "consensus"
)

lrdata = li.mt.inflow(
    adata,
    groupby="fine_cell_type",
    resource=resource,
    use_raw=False
)

print(
    "\nLIANA inflow result shape:"
)

print(
    lrdata.shape
)


###############################################################################
# Visualization: cell types and brain regions
###############################################################################

cell_type_col = "fine_cell_type"
brain_regions = "region"
spatial_key = "spatial"

axs = sc.pl.embedding(
    lrdata,
    basis=spatial_key,
    color=[
        cell_type_col,
        brain_regions
    ],
    s=10,
    ncols=2,
    wspace=0.8,
    show=False
)

for ax in np.ravel(axs):

    ax.invert_yaxis()

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

plt.savefig(
    os.path.join(
        figure_path,
        f"cell_type_brain_region_embedding_{tissue_id}_{bandwidth}.png"
    ),
    dpi=300,
    bbox_inches="tight"
)

plt.close()


###############################################################################
# Example LR interactions
###############################################################################

interactions = [
    "MHb.1^CALM3^TRPC5",
    "MHb.1^NRXN2^NLGN1",
    "MHb.2^NXPH3^NRXN1"
]


###############################################################################
# Plot selected LR interactions
###############################################################################

for interaction in interactions:

    print(
        f"\nPlotting {interaction}..."
    )

    comp = interaction.split("^")

    cell_type = comp[0]
    ligand = comp[1]
    receptor = comp[2]

    safe_interaction = interaction.replace(
        "^",
        "_"
    )


    ###########################################################################
    # 1. Spatial plot of LR interaction score
    ###########################################################################

    ax = sc.pl.embedding(
        lrdata,
        basis=spatial_key,
        color=interaction,
        cmap="viridis_r",
        vmax="p99",
        s=10,
        ncols=2,
        show=False
    )

    ax.invert_yaxis()

    ax.set_aspect(
        "equal",
        adjustable="box"
    )

    plt.savefig(
        os.path.join(
            figure_path,
            f"{tissue_id}_{safe_interaction}_embedding_{bandwidth}.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


    ###########################################################################
    # 2. Spatial plots of ligand and receptor expression
    ###########################################################################

    axs = sc.pl.embedding(
        adata,
        basis=spatial_key,
        color=[
            ligand,
            receptor
        ],
        s=10,
        use_raw=False,
        ncols=2,
        cmap="viridis_r",
        vmax="p99",
        show=False
    )

    for ax in np.ravel(axs):

        ax.invert_yaxis()

        ax.set_aspect(
            "equal",
            adjustable="box"
        )

    plt.savefig(
        os.path.join(
            figure_path,
            f"{tissue_id}_{ligand}_{receptor}_embedding_{bandwidth}.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


    ###########################################################################
    # 3. Violin plot of LR interaction score
    ###########################################################################

    fig, ax = plt.subplots(
        figsize=(14, 5)
    )

    sc.pl.violin(
        lrdata,
        groupby=cell_type_col,
        keys=interaction,
        size=0.5,
        rotation=90,
        ax=ax,
        show=False
    )

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            figure_path,
            f"{tissue_id}_{safe_interaction}_violin_{bandwidth}.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


###############################################################################
# Global interaction summaries
###############################################################################

li.mt.compute_global_specificity(
    lrdata,
    groupby="fine_cell_type",
    use_raw=False,
    verbose=True
)


###############################################################################
# Print top global LR interactions
###############################################################################

print(
    "\nTop 3 global interactions:"
)

print(
    lrdata.uns["global_interactions"]
    .sort_values(
        "lr_mean",
        ascending=False
    )
    .head(3)
)


###############################################################################
# Save LIANA result
###############################################################################

save_path = os.path.join(
    out_path,
    f"lrdata_{tissue_id}_{bandwidth}.h5ad"
)

lrdata.write(
    save_path
)

print(
    f"\nSaved LIANA result:\n{save_path}"
)

print(
    "\nAnalysis completed successfully."
)