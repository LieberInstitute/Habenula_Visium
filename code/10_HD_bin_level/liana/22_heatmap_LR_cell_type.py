#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Heatmap of union(top10 LR pairs across cell types) × cell types
Fill = mean bivariate score (mean_score)

Input:  celltype_specific_interactions_top10.csv
Columns needed: cell_type, interaction, mean_score

Outputs:
- heatmap_top10_union_ctmean.png
- heatmap_top10_union_ctmean_matrix.csv
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist

# -----------------------
# 0) Paths / settings
# -----------------------
infile = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_files/celltype_top10_union_all_celltypes_matrix.csv"
outdir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure"

os.makedirs(outdir, exist_ok=True)

out_png = os.path.join(outdir, "heatmap_top10_union_all_celltypes_from_matrix.png")

# clustering
cluster_rows = True
cluster_cols = True
linkage_method = "average"
distance_metric = "euclidean"

# if NA still exists, choose whether to fill with 0 for plotting
fill_missing_with_zero = False

# figure size
fig_w = 0.6   # width per column
fig_h = 0.25  # height per row
base_w = 6
base_h = 3

# -----------------------
# 1) Load matrix
# -----------------------
mat = pd.read_csv(infile, index_col=0)

# make sure numeric
mat = mat.apply(pd.to_numeric, errors="coerce")

print(f"Matrix shape: {mat.shape}")
print(f"NA cells: {mat.isna().sum().sum()}")

# -----------------------
# 2) Clustering
# -----------------------
# use 0-fill only for clustering
mat_for_cluster = mat.fillna(0.0)

row_order = mat.index
col_order = mat.columns

if cluster_rows and mat.shape[0] > 2:
    Zr = linkage(
        pdist(mat_for_cluster.values, metric=distance_metric),
        method=linkage_method
    )
    row_order = mat.index[leaves_list(Zr)]

if cluster_cols and mat.shape[1] > 2:
    Zc = linkage(
        pdist(mat_for_cluster.values.T, metric=distance_metric),
        method=linkage_method
    )
    col_order = mat.columns[leaves_list(Zc)]

mat_plot = mat.loc[row_order, col_order]

# -----------------------
# 3) Plot heatmap
# -----------------------
if fill_missing_with_zero:
    mat_heatmap = mat_plot.fillna(0.0)
    mask = None
else:
    mat_heatmap = mat_plot
    mask = mat_plot.isna()

sns.set(style="white", context="talk")

n_rows, n_cols = mat_heatmap.shape
plt.figure(figsize=(base_w + fig_w * n_cols, base_h + fig_h * n_rows), dpi=200)

ax = sns.heatmap(
    mat_heatmap,
    mask=mask,
    cmap="OrRd",
    linewidths=0.2,
    linecolor="white",
    cbar_kws={"label": "Mean score"}
)

ax.set_title("Union of top10 LR pairs across cell types")
ax.set_xlabel("Cell type")
ax.set_ylabel("LR pair (interaction)")

plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(out_png, bbox_inches="tight")
plt.close()

print("Done.")
print(f"- Heatmap: {out_png}")



