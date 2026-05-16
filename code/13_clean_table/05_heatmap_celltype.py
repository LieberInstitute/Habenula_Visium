#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from matplotlib.patches import Rectangle

# -----------------------
# 0) Paths / settings
# -----------------------
infile = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_files/celltype_top10_union_all_celltypes_matrix.csv"
outdir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure"

os.makedirs(outdir, exist_ok=True)

out_png = os.path.join(outdir, "heatmap_top10_union_all_celltypes_with_cellwise_shared_unique.png")
out_pdf = os.path.join(outdir, "heatmap_top10_union_all_celltypes_with_cellwise_shared_unique.pdf")
out_ann = os.path.join(outdir, "cellwise_shared_unique_annotation_matrix.csv")

# clustering
cluster_rows = True
cluster_cols = True
linkage_method = "average"
distance_metric = "euclidean"

# if NA still exists, choose whether to fill with 0 for plotting
fill_missing_with_zero = False

# figure size
fig_w = 0.6
fig_h = 0.25
base_w = 6
base_h = 3

# -----------------------
# 0b) Cell-wise annotation settings
# -----------------------
eps = 1e-8

# absolute activity cutoff on original scale
abs_cutoff = 0.05

# relative activity cutoff within each row
rel_cutoff = 0.75

# uniqueness threshold:
# top value must be >= unique_ratio * second-highest value
unique_ratio = 1.3

# colors
unique_edgecolor = "#E45756"  # red
shared_edgecolor = "#4C78A8"  # blue

# -----------------------
# 1) Load matrix
# -----------------------
mat = pd.read_csv(infile, index_col=0)
mat = mat.apply(pd.to_numeric, errors="coerce")

print(f"Matrix shape: {mat.shape}")
print(f"NA cells: {mat.isna().sum().sum()}")

# -----------------------
# 2) Clustering (same as before)
# -----------------------
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
# 3) Cell-wise shared / unique annotation
# -----------------------
mat_nonneg = mat.fillna(0.0).clip(lower=0)

annotation = pd.DataFrame(
    "other",
    index=mat_nonneg.index,
    columns=mat_nonneg.columns
)

for lr in mat_nonneg.index:
    x = mat_nonneg.loc[lr].values.astype(float)

    if np.all(x <= 0):
        continue

    row_max = np.max(x)
    if row_max <= 0:
        continue

    # relative-to-row-max scores
    r = x / (row_max + eps)

    # active cell types for this LR
    active = (x >= abs_cutoff) & (r >= rel_cutoff)
    active_idx = np.where(active)[0]
    active_n = len(active_idx)

    if active_n == 0:
        continue

    # sort to get top1 / top2
    sort_idx = np.argsort(x)[::-1]
    top1_idx = sort_idx[0]
    top1_val = x[top1_idx]
    top2_val = x[sort_idx[1]] if len(sort_idx) > 1 else 0.0

    # classify
    if active_n == 1:
        annotation.iloc[mat_nonneg.index.get_loc(lr), active_idx[0]] = "unique"
    else:
        # if top cell type strongly dominates, mark it unique and the rest shared
        if top2_val <= 0:
            annotation.iloc[mat_nonneg.index.get_loc(lr), top1_idx] = "unique"
            for j in active_idx:
                if j != top1_idx:
                    annotation.iloc[mat_nonneg.index.get_loc(lr), j] = "shared"
        elif (top1_val / (top2_val + eps)) >= unique_ratio:
            annotation.iloc[mat_nonneg.index.get_loc(lr), top1_idx] = "unique"
            for j in active_idx:
                if j != top1_idx:
                    annotation.iloc[mat_nonneg.index.get_loc(lr), j] = "shared"
        else:
            # several cell types similarly high -> shared
            for j in active_idx:
                annotation.iloc[mat_nonneg.index.get_loc(lr), j] = "shared"

# save full annotation matrix
annotation.to_csv(out_ann)

print("\nCell-wise annotation counts:")
print(annotation.stack().value_counts())
print(f"Saved annotation matrix: {out_ann}")

# reorder annotation to match plotted heatmap
annotation_plot = annotation.loc[row_order, col_order]

# -----------------------
# 4) Plot heatmap
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

# -----------------------
# 5) Overlay cell-wise unique/shared markers
# -----------------------
for i, lr in enumerate(annotation_plot.index):
    for j, ct in enumerate(annotation_plot.columns):
        label = annotation_plot.loc[lr, ct]

        if label == "unique":
            rect = Rectangle(
                (j, i), 1, 1,
                fill=False,
                edgecolor=unique_edgecolor,
                linewidth=1.5
            )
            ax.add_patch(rect)

        elif label == "shared":
            rect = Rectangle(
                (j, i), 1, 1,
                fill=False,
                edgecolor=shared_edgecolor,
                linewidth=1.0
            )
            ax.add_patch(rect)

# legend
handles = [
    Rectangle((0, 0), 1, 1, fill=False, edgecolor=unique_edgecolor, linewidth=1.5, label="unique"),
    Rectangle((0, 0), 1, 1, fill=False, edgecolor=shared_edgecolor, linewidth=1.0, label="shared")
]

ax.legend(
    handles=handles,
    title="Cell-wise label",
    bbox_to_anchor=(1.25, 1),
    loc="upper left",
    borderaxespad=0
)

plt.subplots_adjust(right=0.82)
plt.tight_layout()
plt.savefig(out_png, bbox_inches="tight")
plt.savefig(out_pdf, bbox_inches="tight")
plt.close()


print("\nDone.")
print(f"- Heatmap: {out_png}")
print(f"- Heatmap PDF: {out_pdf}")
print(f"- Annotation matrix: {out_ann}")