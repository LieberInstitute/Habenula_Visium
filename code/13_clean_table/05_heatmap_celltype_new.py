#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

out_png = os.path.join(outdir, "heatmap_top10_union_all_celltypes_markX_rowmean_colcluster.png")
out_pdf = os.path.join(outdir, "heatmap_top10_union_all_celltypes_markX_rowmean_colcluster.pdf")

# -----------------------
# 1) Settings
# -----------------------
fill_missing_with_zero = False
top_n = 10

# column clustering only
cluster_cols = True
linkage_method = "average"
distance_metric = "euclidean"

# figure size
fig_w = 0.6
fig_h = 0.25
base_w = 6
base_h = 3

# X marker style
x_color = "black"
x_fontsize = 9
x_fontweight = "bold"

# heatmap style
cmap = "OrRd"

# -----------------------
# 2) Load matrix
# -----------------------
mat = pd.read_csv(infile, index_col=0)
mat = mat.apply(pd.to_numeric, errors="coerce")

print(f"Matrix shape: {mat.shape}")
print(f"NA cells: {mat.isna().sum().sum()}")

# -----------------------
# 3) Row order: sort by mean descending
# -----------------------
row_means = mat.mean(axis=1, skipna=True)
row_order = row_means.sort_values(ascending=False).index

# -----------------------
# 4) Column order: clustering
# -----------------------
mat_for_cluster = mat.fillna(0.0)

col_order = mat.columns
if cluster_cols and mat.shape[1] > 2:
    Zc = linkage(
        pdist(mat_for_cluster.values.T, metric=distance_metric),
        method=linkage_method
    )
    col_order = mat.columns[leaves_list(Zc)]

# reordered matrix
mat_plot = mat.loc[row_order, col_order]

# -----------------------
# 5) Mark top 10 within each column with X
# -----------------------
marker = pd.DataFrame(
    "",
    index=mat_plot.index,
    columns=mat_plot.columns
)

for ct in mat_plot.columns:
    col = mat_plot[ct].dropna()

    if col.empty:
        continue

    top_idx = col.sort_values(ascending=False).head(top_n).index
    marker.loc[top_idx, ct] = "X"

# -----------------------
# 6) Prepare heatmap matrix
# -----------------------
if fill_missing_with_zero:
    mat_heatmap = mat_plot.fillna(0.0)
    mask = None
else:
    mat_heatmap = mat_plot
    mask = mat_plot.isna()

# -----------------------
# 7) Plot heatmap
# -----------------------
sns.set(style="white", context="talk")

n_rows, n_cols = mat_heatmap.shape
plt.figure(figsize=(base_w + fig_w * n_cols, base_h + fig_h * n_rows), dpi=200)

ax = sns.heatmap(
    mat_heatmap,
    mask=mask,
    cmap=cmap,
    linewidths=0.2,
    linecolor="white",
    cbar_kws={"label": "Mean score"}
)

ax.set_title("Union of top10 LR pairs across cell types")
ax.set_xlabel("Cell type")
ax.set_ylabel("LR pair (interaction)")

plt.xticks(rotation=90)
plt.yticks(rotation=0)

# -----------------------
# 8) Overlay X on top10 cells
# -----------------------
for i, lr in enumerate(marker.index):
    for j, ct in enumerate(marker.columns):
        if marker.loc[lr, ct] == "X":
            ax.text(
                j + 0.5,
                i + 0.5,
                "X",
                ha="center",
                va="center",
                color=x_color,
                fontsize=x_fontsize,
                fontweight=x_fontweight
            )

plt.tight_layout()
plt.savefig(out_png, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()

print("\nDone.")
print(f"- Heatmap PNG: {out_png}")
print(f"- Heatmap PDF: {out_pdf}")