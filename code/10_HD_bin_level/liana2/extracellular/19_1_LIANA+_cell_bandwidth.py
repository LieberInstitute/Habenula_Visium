# Try difference bandwith to see the results

import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
import muon as mu
import anndata as ad

import seaborn as sns
import matplotlib.pyplot as plt
import plotnine as p9

import squidpy as sq
import os
from pyhere import here
import decoupler as dc
import plotnine as p9

import matplotlib.pyplot as plt
import json

# Choose bandwidth
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--bandwidth", type=float, required=True)
args = parser.parse_args()

bandwidth = args.bandwidth

print("bandwidth =", bandwidth)

in_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana2','extracellular'
)
os.makedirs(os.path.join(out_path), exist_ok=True)
figure_path = here("plots", "10_HD_bin_level", "no_secondary", "liana2","extracellular")
os.makedirs(figure_path, exist_ok=True)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
sample_info = pd.read_csv(sample_info_path)
sample_id = sample_info['sample_id'].iloc[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]
tissue_id = sample_info['tissue_id'].iloc[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

scale_json = here(
    'processed-data', '01_spaceranger', "five_samples_10_2025", sample_id,
    'outs', 'binned_outputs', 'square_008um', 'spatial',
    'scalefactors_json.json'
)

file_path = os.path.join(in_path, f"adata/adata_extracellular_withcelltype_newcluster.h5ad")
adata = sc.read(file_path)
adata.var.index.name = None

###### run separately for each sample ######
adata.obs = adata.obs.rename(columns={"tissue_section": "tissue_id"})
adata = adata[adata.obs["tissue_id"].isin([tissue_id])].copy()

adata = adata[:, ~adata.var_names.str.startswith("DEPRECATED_")].copy()
adata.var_names = adata.var_names.astype(str)
adata.var_names_make_unique()

adata = adata[adata.obs["fine_cell_type"] != "Drop"].copy()

# figure
ax = sc.pl.embedding(
    adata,
    basis="spatial",
    color=["fine_cell_type"],
    wspace=0.4,
    s=5,
    show=False
)

ax.invert_yaxis()
ax.set_aspect("equal", adjustable="box")

plt.savefig(
    os.path.join(figure_path, f"{tissue_id}_fine_cell_type_spatial.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Basic Prep and QC
## filter cells and genes
sc.pp.filter_cells(adata, min_genes=10)
sc.pp.filter_genes(adata, min_cells=3) 

adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# calculate spatial connectivity matrix Wij
plot, df = li.ut.query_bandwidth(
    coordinates=adata.obsm["spatial"],
    start=0,
    end=adata.n_obs-1,
    interval_n=40
)

ymin = int(df.neighbours.min())
ymax = int(df.neighbours.max())
step = max(1, round((ymax - ymin) / 6))

p = (
    plot
    + p9.geom_vline(xintercept=bandwidth, linetype="dashed", color="red")
    + p9.annotate(
        "text",
        x=bandwidth,
        y=ymax,
        label=f"chosen = {bandwidth:.0f} µm",
        angle=90,
        va="bottom",
        ha="right",
        color="red"
    )
    + p9.scale_y_continuous(
        breaks=list(range(ymin, ymax + 1, step))
    )
)

p.save(
    filename=os.path.join(
        figure_path,
        f"{tissue_id}_spatial_connectivity_bandwidth_{bandwidth}.png"
    ),
    dpi=300,
    width=6,
    height=4,
    units="in"
)

li.ut.spatial_neighbors(adata=adata, bandwidth=bandwidth, spatial_key="spatial", max_neighbours=adata.n_obs-1)

idx = min(5500, adata.n_obs - 1)
# Visualize proximity graph

p = (
    li.pl.connectivity(
    adata,
    idx=idx,
    size=0.01,
    # figure_size=(6, 10),
    spatial_key="spatial",
    return_fig=True)
     + p9.coord_fixed()
)

p.save(
    filename=os.path.join(
        figure_path,
        f"{tissue_id}_connectivity_idx{idx}_{bandwidth}.png"
    ),
    dpi=300,
    width=6,
    height=10,
    units="in"
)

## Find spatially variable genes (SVGs) using Morans’I

# sq.gr.spatial_autocorr(adata, mode='moran', use_raw=False, show_progress_bar=True)

## Filter by spatially variable genes (SVGs)
## Check how many spatially variable genes (SVGs)
# svgs = adata.uns['moranI'].index[(adata.uns['moranI']['pval_norm_fdr_bh'] < 0.05) & (adata.uns['moranI']['I'] > 0.01)]
# len(svgs)

# adata = adata[:, svgs]

# Instead of using this,we will use our SVG list
# list_path = here('processed-data', '10_HD_bin_level', 'no_secondary', 'nnSVG_out', 'merged_SVGs.txt')
# with open(list_path, 'r') as f:
#     svgs = f.read().splitlines()

# svg_mask = adata.var["gene_id"].isin(svgs)
# adata = adata[:, svg_mask]

# Compute inflow score
resource = li.rs.select_resource("consensus")
lrdata = li.mt.inflow(adata,
                      groupby='fine_cell_type',
                      resource=resource,
                      use_raw=False)
                                  
lrdata.shape

# visualization
## Define variables
cell_type_col = "fine_cell_type"
brain_regions = "region"
spatial_key = "spatial"

axs = sc.pl.embedding(
    lrdata,
    basis=spatial_key,
    color=[cell_type_col, brain_regions],
    s=10,
    ncols=2,
    wspace=0.8,
    show=False
)

for ax in np.ravel(axs):
    ax.invert_yaxis()
    ax.set_aspect("equal", adjustable="box")

plt.savefig(
    os.path.join(figure_path, f"cell_type_brain_region_embedding_{tissue_id}_{bandwidth}.png"),
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# interaction = 'MHb.1^ITGAV^THY1'
# interaction = 'MHb.1^APP^GPC1'
# interaction = 'MHb.1^LIN7C^HTR2C'

interactions = [
    "MHb.1^CALM3^TRPC5",
    "MHb.1^NRXN2^NLGN1",
    "MHb.2^NXPH3^NRXN1"
]

for interaction in interactions:

    print(f"Plotting {interaction}...")

    comp = interaction.split("^")
    cell_type = comp[0]
    ligand = comp[1]
    receptor = comp[2]

    safe_interaction = interaction.replace("^", "_")

    # -----------------------------
    # 1. Plot LR interaction score
    # -----------------------------
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
    ax.set_aspect("equal", adjustable="box")

    plt.savefig(
        os.path.join(
            figure_path,
            f"{tissue_id}_{safe_interaction}_embedding_{bandwidth}.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # -----------------------------
    # 2. Plot ligand and receptor expression
    # -----------------------------
    axs = sc.pl.embedding(
        adata,
        basis=spatial_key,
        color=[ligand, receptor],
        s=10,
        use_raw=False,
        ncols=2,
        cmap="viridis_r",
        vmax="p99",
        show=False
    )

    for ax in np.ravel(axs):
        ax.invert_yaxis()
        ax.set_aspect("equal", adjustable="box")

    plt.savefig(
        os.path.join(
            figure_path,
            f"{tissue_id}_{ligand}_{receptor}_embedding_{bandwidth}.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

    # -----------------------------
    # 3. Violin plot of LR interaction score
    # -----------------------------
    fig, ax = plt.subplots(figsize=(14, 5))

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

# Global Summaries

li.mt.compute_global_specificity(lrdata, groupby='fine_cell_type', use_raw=False, verbose=True)
lrdata.uns['global_interactions'].sort_values("lr_mean", ascending=False).head(3)

save_path = os.path.join(out_path, f"lrdata_{tissue_id}_{bandwidth}.h5ad")
lrdata.write(save_path)