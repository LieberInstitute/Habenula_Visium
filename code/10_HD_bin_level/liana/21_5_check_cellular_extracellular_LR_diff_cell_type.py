#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Spearman cellular vs extracellular LR rank plots by cell type
+ Outlier detection (distance from trend line within each cell type)
+ Union-of-outliers heatmap (distance + optional binary)

Input:
- celltype_specific_interactions_all.csv (cellular)
- celltype_specific_interactions_all.csv (extracellular)

Output:
- <plots_dir>/rank_rank_plots/rank_rank_<cell_type>.png
- <plots_dir>/heatmap_outlier_union_distance.png
- <plots_dir>/heatmap_outlier_union_binary.png (optional but produced)
- <plots_dir>/outlier_union_table.csv (optional: outlier table)
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


# -----------------------
# 0) Config
# -----------------------
plots_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure"

cell_path = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_specific_interactions/celltype_specific_interactions_all.csv"
extra_path = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_specific_interactions_extracellular/celltype_specific_interactions_all.csv"

value_col = "ct_mean"
rank_method = "average"   # "average" or "dense" if many ties
min_points_per_celltype = 10
z_thresh = 3.0            # robust z threshold for outliers
top_annot = 8             # number of outliers to annotate per plot

# If you want deterministic plot style
plt.rcParams.update({"figure.max_open_warning": 0})

os.makedirs(plots_dir, exist_ok=True)


# -----------------------
# 1) Helpers
# -----------------------
def prep_rank(df, value_col="ct_mean", dataset_name="cell", rank_method="average"):
    """
    Deduplicate/aggregate by (cell_type, interaction),
    then compute within-cell_type rank by descending value_col.
    """
    df2 = (df
           .groupby(["cell_type", "interaction"], as_index=False)[value_col]
           .mean()
           .rename(columns={value_col: f"{value_col}_{dataset_name}"}))

    df2[f"rank_{dataset_name}"] = (
        df2.groupby("cell_type")[f"{value_col}_{dataset_name}"]
           .rank(ascending=False, method=rank_method)
    )
    return df2


def fit_line_and_distance(x, y):
    """
    Fit y = a + b x and compute perpendicular distance to line:
    |b x - y + a| / sqrt(b^2 + 1)

    Uses Theil–Sen robust regression if sklearn is available,
    otherwise falls back to OLS polyfit.
    """
    x = np.asarray(x).reshape(-1)
    y = np.asarray(y).reshape(-1)

    try:
        from sklearn.linear_model import TheilSenRegressor
        model = TheilSenRegressor(random_state=0)
        model.fit(x.reshape(-1, 1), y)
        b = float(model.coef_[0])
        a = float(model.intercept_)
    except Exception:
        b, a = np.polyfit(x, y, 1)

    dist = np.abs(b * x - y + a) / np.sqrt(b**2 + 1.0)
    return a, b, dist


def robust_z(x):
    """
    Robust z-score using MAD:
    (x - median) / (1.4826 * MAD)
    """
    x = np.asarray(x)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    denom = 1.4826 * mad + 1e-12
    return (x - med) / denom


def detect_outliers_one_celltype(df_sub, z_thresh=3.0):
    """
    For a single cell_type subset (merged), compute:
    - fitted line params
    - distance to line
    - robust z of distance
    - outlier flag (robust z > z_thresh)
    """
    a, b, dist = fit_line_and_distance(df_sub["rank_cell"], df_sub["rank_extra"])
    rz = robust_z(dist)

    out = df_sub.copy()
    out["line_a"] = a
    out["line_b"] = b
    out["dist"] = dist
    out["dist_rz"] = rz
    out["outlier"] = (rz > z_thresh)
    return out


def plot_rank_rank_one_celltype(df_sub, outdir, top_annot=8):
    """
    Scatter rank_cell vs rank_extra, overlay trend line,
    color outliers, annotate top outliers by distance.
    """
    os.makedirs(outdir, exist_ok=True)
    ct = df_sub["cell_type"].iloc[0]

    rho, p = spearmanr(df_sub["rank_cell"], df_sub["rank_extra"])
    a = df_sub["line_a"].iloc[0]
    b = df_sub["line_b"].iloc[0]

    x = df_sub["rank_cell"].to_numpy()
    y = df_sub["rank_extra"].to_numpy()

    xs = np.linspace(x.min(), x.max(), 200)
    ys = a + b * xs

    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)

    inl = df_sub[~df_sub["outlier"]]
    out = df_sub[df_sub["outlier"]]

    ax.scatter(inl["rank_cell"], inl["rank_extra"], s=10, alpha=0.6)
    ax.scatter(out["rank_cell"], out["rank_extra"], s=18, alpha=0.9)
    ax.plot(xs, ys, linewidth=1)

    ax.set_title(f"{ct} | Spearman rho={rho:.3f} (p={p:.1e}) | n={len(df_sub)}")
    ax.set_xlabel("Cellular rank (ct_mean, within cell type)")
    ax.set_ylabel("Extracellular rank (ct_mean, within cell type)")

    if out.shape[0] > 0:
        out2 = out.sort_values("dist", ascending=False).head(top_annot)
        for _, r in out2.iterrows():
            ax.text(r["rank_cell"], r["rank_extra"], str(r["interaction"]), fontsize=6)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"rank_rank_{ct}.png"))
    plt.close(fig)


def cluster_columns(mat, method="average", metric="euclidean"):
    """
    Cluster columns (cell types) and return matrix with reordered columns.
    """
    if mat.shape[1] <= 2:
        return mat
    Z = linkage(pdist(mat.T, metric=metric), method=method)
    order = leaves_list(Z)
    return mat.iloc[:, order]


def save_heatmap(mat, outfile, title, cbar_label, figsize_scale=(0.35, 0.18)):
    """
    Simple matplotlib heatmap (no seaborn dependency).
    """
    h, w = mat.shape
    fig_w = figsize_scale[0] * w + 4
    fig_h = figsize_scale[1] * h + 4

    plt.figure(figsize=(fig_w, fig_h), dpi=150)
    plt.imshow(mat.values, aspect="auto")
    plt.yticks(np.arange(h), mat.index, fontsize=6)
    plt.xticks(np.arange(w), mat.columns, rotation=90, fontsize=7)
    plt.colorbar(label=cbar_label)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


# -----------------------
# 2) Load
# -----------------------
df_cell = pd.read_csv(cell_path)
df_extra = pd.read_csv(extra_path)

required_cols = {"cell_type", "interaction", value_col}
if not required_cols.issubset(df_cell.columns) or not required_cols.issubset(df_extra.columns):
    raise ValueError(f"Input CSVs must contain columns: {required_cols}")

# -----------------------
# 3) Rank within cell_type (by ct_mean)
# -----------------------
cell_r = prep_rank(df_cell, value_col=value_col, dataset_name="cell", rank_method=rank_method)
extra_r = prep_rank(df_extra, value_col=value_col, dataset_name="extra", rank_method=rank_method)

df_join = cell_r.merge(extra_r, on=["cell_type", "interaction"], how="inner")

# -----------------------
# 4) Outlier detection within each cell_type
# -----------------------
all_res = []
for ct, sub in df_join.groupby("cell_type"):
    if sub.shape[0] < min_points_per_celltype:
        continue
    all_res.append(detect_outliers_one_celltype(sub, z_thresh=z_thresh))

if len(all_res) == 0:
    raise RuntimeError("No cell types passed the min_points_per_celltype filter.")

df_res = pd.concat(all_res, ignore_index=True)

# Save full merged table with outlier info
df_res_outfile = os.path.join(plots_dir, "cellular_vs_extracellular_rank_with_outliers.csv")
df_res.to_csv(df_res_outfile, index=False)

# -----------------------
# 5) Rank-rank plots per cell type
# -----------------------
rank_plot_dir = os.path.join(plots_dir, "rank_rank_plots")
for ct, sub in df_res.groupby("cell_type"):
    plot_rank_rank_one_celltype(sub, outdir=rank_plot_dir, top_annot=top_annot)

# -----------------------
# 6) Union-of-outliers heatmaps (distance + binary)
# -----------------------
union_interactions = (df_res.loc[df_res["outlier"], "interaction"]
                      .drop_duplicates()
                      .sort_values()
                      .tolist())

if len(union_interactions) == 0:
    # still make an empty note file
    with open(os.path.join(plots_dir, "NO_OUTLIERS_FOUND.txt"), "w") as f:
        f.write(f"No outliers found with z_thresh={z_thresh}\n")
else:
    df_u = df_res[df_res["interaction"].isin(union_interactions)].copy()

    # Distance matrix: fill non-present with 0 (or use NaN if you prefer)
    mat_dist = df_u.pivot_table(index="interaction", columns="cell_type", values="dist", aggfunc="max")
    mat_dist = mat_dist.fillna(0.0)

    # Binary matrix (avoid pandas FutureWarning by making 0/1 before pivot)
    df_u["outlier01"] = df_u["outlier"].astype("int8")
    mat_bin = df_u.pivot_table(index="interaction", columns="cell_type", values="outlier01", aggfunc="max")
    mat_bin = mat_bin.fillna(0).astype("int8")

    # Cluster columns (cell types)
    mat_dist_c = cluster_columns(mat_dist)
    mat_bin_c = mat_bin[mat_dist_c.columns]  # keep same column order as dist

    # Save heatmaps
    save_heatmap(
        mat_dist_c,
        outfile=os.path.join(plots_dir, "heatmap_outlier_union_distance.png"),
        title="Outlier union heatmap (distance from trend line)",
        cbar_label="Distance to trend line"
    )

    save_heatmap(
        mat_bin_c,
        outfile=os.path.join(plots_dir, "heatmap_outlier_union_binary.png"),
        title="Outlier union heatmap (binary)",
        cbar_label="Outlier (0/1)"
    )

    # Optional: save the union matrix tables
    mat_dist_c.to_csv(os.path.join(plots_dir, "heatmap_outlier_union_distance.csv"))
    mat_bin_c.to_csv(os.path.join(plots_dir, "heatmap_outlier_union_binary.csv"))

    # Optional: save union outlier rows (long form)
    df_u.sort_values(["cell_type", "dist"], ascending=[True, False]) \
        .to_csv(os.path.join(plots_dir, "outlier_union_table.csv"), index=False)


print("Done.")
print(f"- Rank-rank plots: {rank_plot_dir}")
print(f"- Merged outlier table: {df_res_outfile}")
print(f"- Heatmaps saved under: {plots_dir}")
