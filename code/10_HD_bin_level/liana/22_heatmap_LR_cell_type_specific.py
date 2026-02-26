#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Heatmap of union(top10 LR pairs across cell types) × cell types
Fill = mean bivariate score (abs_z)

Input:  celltype_specific_interactions_top10.csv
Columns needed: cell_type, interaction, abs_z

Outputs:
- heatmap_top10_union_ctmean.png
- heatmap_top10_union_ctmean_matrix.csv
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


# -----------------------
# 0) Paths / settings
# -----------------------
infile = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_specific_interactions/celltype_specific_interactions_top10.csv" 
outdir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure"

out_png = os.path.join(outdir, "heatmap_top10_union_ctmean_specific.png")
out_csv = os.path.join(outdir, "heatmap_top10_union_ctmean_matrix_specific.csv")

cluster_rows = True
cluster_cols = True
linkage_method = "average"
distance_metric = "euclidean"

# Heatmap appearance
fig_w = 0.6   # width per column
fig_h = 0.25  # height per row
base_w = 6
base_h = 3


# -----------------------
# 1) Load
# -----------------------
os.makedirs(outdir, exist_ok=True)
df = pd.read_csv(infile)

need = {"cell_type", "interaction", "abs_z"}
missing = need - set(df.columns)
if missing:
    raise ValueError(f"Missing columns in {infile}: {missing}")

df = df[["cell_type", "interaction", "abs_z"]].copy()

# If there are duplicates, keep the max abs_z (or change to mean if you prefer)
df = (df.groupby(["cell_type", "interaction"], as_index=False)["abs_z"]
        .max())

# -----------------------
# 2) Union of top10 across cell types -> matrix
# -----------------------
mat = df.pivot_table(
    index="interaction",
    columns="cell_type",
    values="abs_z",
    aggfunc="max"
)

# Save raw matrix (NaN means that LR pair is not in top10 of that cell type)
mat.to_csv(out_csv)

# -----------------------
# 3) Optional clustering (rows/cols)
#    (Use 0-fill only for clustering distance calc; keep NaN for plotting mask)
# -----------------------
mat_for_cluster = mat.fillna(0.0)

row_order = mat.index
col_order = mat.columns

if cluster_rows and mat.shape[0] > 2:
    Zr = linkage(pdist(mat_for_cluster.values, metric=distance_metric), method=linkage_method)
    row_order = mat.index[leaves_list(Zr)]

if cluster_cols and mat.shape[1] > 2:
    Zc = linkage(pdist(mat_for_cluster.values.T, metric=distance_metric), method=linkage_method)
    col_order = mat.columns[leaves_list(Zc)]

mat_plot = mat.loc[row_order, col_order]

# -----------------------
# 4) Plot heatmap
# -----------------------
sns.set(style="white", context="talk")

n_rows, n_cols = mat_plot.shape
plt.figure(figsize=(base_w + fig_w * n_cols, base_h + fig_h * n_rows), dpi=200)

# Mask NaN so missing pairs are blank
mask = mat_plot.isna()

ax = sns.heatmap(
    mat_plot,
    mask=mask,
    linewidths=0.2,
    linecolor="white",
    cbar_kws={"label": "Mean bivariate score (abs_z)"}
)

ax.set_title("Union of top10 LR pairs across cell types (fill = abs_z)")
ax.set_xlabel("Cell type specifc")
ax.set_ylabel("LR pair (interaction)")

plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig(out_png)
plt.close()

print("Done.")
print(f"- Heatmap: {out_png}")
print(f"- Matrix CSV: {out_csv}")








