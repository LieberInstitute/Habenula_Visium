"""
Rank–Rank comparison between cellular vs extracellular interaction ranks,
outlier detection per cell type, and two-axis clustered heatmaps.

Inputs:
- celltype_specific_interactions_all.csv (cellular)
- celltype_specific_interactions_all.csv (extracellular)

Outputs:
- <plots_dir>/rank_rank_plots/rank_rank_<cell_type>.png
- <plots_dir>/cellular_vs_extracellular_rank_with_outliers.csv
- <plots_dir>/heatmap_outlier_union_distance.png
- <plots_dir>/heatmap_outlier_union_binary.png
- <plots_dir>/heatmap_outlier_union_signed_distance.png
- <plots_dir>/heatmap_outlier_union_signed_binary.png
- <plots_dir>/heatmap_outlier_union_distance.csv
- <plots_dir>/heatmap_outlier_union_binary.csv
- <plots_dir>/heatmap_outlier_union_signed_distance.csv
- <plots_dir>/heatmap_outlier_union_signed_binary.csv
- <plots_dir>/outlier_union_table.csv
- <plots_dir>/NO_OUTLIERS_FOUND.txt (if no outliers)
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
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
rank_method = "average"   # {"average", "dense"}
min_points_per_celltype = 10
z_thresh = 3.5            # robust z threshold for outliers (based on ABS distance)
top_annot = 8             # number of outliers to annotate per scatter plot

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
    df2 = (
        df.groupby(["cell_type", "interaction"], as_index=False)[value_col]
          .mean()
          .rename(columns={value_col: f"{value_col}_{dataset_name}"})
    )

    df2[f"rank_{dataset_name}"] = (
        df2.groupby("cell_type")[f"{value_col}_{dataset_name}"]
           .rank(ascending=False, method=rank_method)
    )
    return df2


def fit_line_and_distance(x, y):
    """
    Fit y = a + b x and compute perpendicular ABS distance to the line:
    |b x - y + a| / sqrt(b^2 + 1)

    Uses Theil–Sen robust regression if sklearn is available,
    otherwise falls back to NumPy OLS polyfit.
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

    dist = np.abs(b * x - y + a) / np.sqrt(b**2 + 1.0)  # always >= 0
    return a, b, dist


def robust_z(x):
    """
    Robust z-score using MAD:
    (x - median) / (1.4826 * MAD)
    If MAD == 0, return zeros (avoid exploding z).
    """
    x = np.asarray(x, dtype=float)
    med = np.median(x)
    mad = np.median(np.abs(x - med))
    if mad < 1e-12:
        return np.zeros_like(x, dtype=float)
    return (x - med) / (1.4826 * mad)


def detect_outliers_one_celltype(df_sub, z_thresh=3.0):
    """
    For a single cell_type subset (merged), compute:
    - fitted line params (a, b)
    - abs distance to the fitted line (dist): used for outlier detection
    - signed residual and signed distance: used for direction visualization
    - robust z of dist
    - outlier flag (robust z > z_thresh)
    """
    a, b, dist = fit_line_and_distance(df_sub["rank_cell"], df_sub["rank_extra"])

    x = df_sub["rank_cell"].to_numpy()
    y = df_sub["rank_extra"].to_numpy()
    yhat = a + b * x
    residual = y - yhat
    dist_signed = residual / np.sqrt(b**2 + 1.0)

    rz = robust_z(dist)

    out = df_sub.copy()
    out["line_a"] = a
    out["line_b"] = b
    out["dist"] = dist
    out["residual"] = residual
    out["dist_signed"] = dist_signed
    out["direction"] = np.where(residual < 0, "extracellular_leaning", "cellular_leaning")
    out["dist_rz"] = rz
    out["outlier"] = (rz > z_thresh)
    return out


def plot_rank_rank_one_celltype(df_sub, outdir, top_annot=8):
    """
    Scatter rank_cell vs rank_extra, overlay trend line,
    color by signed distance (neg=blue, pos=orange-ish),
    and draw a black ring for outliers.
    """
    os.makedirs(outdir, exist_ok=True)
    ct = df_sub["cell_type"].iloc[0]

    rho, p = spearmanr(df_sub["rank_cell"], df_sub["rank_extra"])
    a = df_sub["line_a"].iloc[0]
    b = df_sub["line_b"].iloc[0]

    x = df_sub["rank_cell"].to_numpy()
    y = df_sub["rank_extra"].to_numpy()
    sd = df_sub["dist_signed"].to_numpy()

    xs = np.linspace(x.min(), x.max(), 200)
    ys = a + b * xs

    vmax = np.nanquantile(np.abs(sd), 0.99)
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1.0
    norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)

    fig, ax = plt.subplots(figsize=(6, 5), dpi=150)

    sc = ax.scatter(
        x, y,
        c=sd,
        cmap="coolwarm",
        norm=norm,
        s=10,
        alpha=0.75,
        linewidths=0
    )

    out = df_sub[df_sub["outlier"]]
    if out.shape[0] > 0:
        ax.scatter(
            out["rank_cell"], out["rank_extra"],
            s=45,
            facecolors="none",
            edgecolors="black",
            linewidths=0.9,
            alpha=0.95,
            label="Outliers"
        )
        out2 = out.assign(absdist=np.abs(out["dist"])).sort_values("absdist", ascending=False).head(top_annot)
        for _, r in out2.iterrows():
            ax.text(r["rank_cell"], r["rank_extra"], str(r["interaction"]), fontsize=6)

    ax.plot(xs, ys, linewidth=1, color="black", alpha=0.7, label="Trend")

    ax.set_title(f"{ct} | Spearman rho={rho:.3f} (p={p:.1e}) | n={len(df_sub)}")
    ax.set_xlabel("Cellular rank (ct_mean, within cell type)")
    ax.set_ylabel("Extracellular rank (ct_mean, within cell type)")

    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Signed distance to trend line\n(neg=extracellular-leaning, pos=cellular-leaning)")

    if out.shape[0] > 0:
        ax.legend(frameon=False, fontsize=7)

    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"rank_rank_{ct}.png"))
    plt.close(fig)


def _safe_pdist(X, metric):
    d = pdist(X, metric=metric)
    if np.isnan(d).any():
        warnings.warn(f"NaNs encountered in pdist(metric='{metric}'); falling back to 'euclidean'.")
        d = pdist(X, metric='euclidean')
    return d


def cluster_rows_and_cols(mat, method="average", metric="correlation"):
    r_order = np.arange(mat.shape[0])
    c_order = np.arange(mat.shape[1])

    if mat.shape[0] > 2:
        Zr = linkage(_safe_pdist(mat.values, metric=metric), method=method)
        r_order = leaves_list(Zr)

    if mat.shape[1] > 2:
        Zc = linkage(_safe_pdist(mat.values.T, metric=metric), method=method)
        c_order = leaves_list(Zc)

    mat2 = mat.iloc[r_order, c_order]
    return mat2, r_order, c_order


def save_heatmap(mat, outfile, title, cbar_label, cmap=None, center0=False, figsize_scale=(0.35, 0.18)):
    h, w = mat.shape
    fig_w = figsize_scale[0] * w + 4
    fig_h = figsize_scale[1] * h + 4

    plt.figure(figsize=(fig_w, fig_h), dpi=150)

    if center0:
        vmax = np.nanquantile(np.abs(mat.values), 0.99)
        if not np.isfinite(vmax) or vmax <= 0:
            vmax = 1.0
        norm = TwoSlopeNorm(vcenter=0.0, vmin=-vmax, vmax=vmax)
    else:
        norm = None

    plt.imshow(mat.values, aspect="auto", cmap=cmap, norm=norm)
    plt.yticks(np.arange(h), mat.index, fontsize=6)
    plt.xticks(np.arange(w), mat.columns, rotation=90, fontsize=7)
    plt.colorbar(label=cbar_label)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()


def max_abs_agg(s):
    """Aggregator: return the value with the largest absolute magnitude (keeps negative)."""
    s = pd.Series(s).dropna()
    if s.empty:
        return np.nan
    i = np.argmax(np.abs(s.to_numpy()))
    return s.iloc[i]

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

# Save full merged table with outlier + sign info
df_res_outfile = os.path.join(plots_dir, "cellular_vs_extracellular_rank_with_outliers.csv")
df_res.to_csv(df_res_outfile, index=False)

# -----------------------
# 5) Rank–rank plots per cell type (colored by sign)
# -----------------------
rank_plot_dir = os.path.join(plots_dir, "rank_rank_plots")
for ct, sub in df_res.groupby("cell_type"):
    plot_rank_rank_one_celltype(sub, outdir=rank_plot_dir, top_annot=top_annot)

# -----------------------
# 6) Union-of-outliers heatmaps
# -----------------------
union_interactions = (
    df_res.loc[df_res["outlier"], "interaction"]
          .drop_duplicates()
          .sort_values()
          .tolist()
)

if len(union_interactions) == 0:
    with open(os.path.join(plots_dir, "NO_OUTLIERS_FOUND.txt"), "w") as f:
        f.write(f"No outliers found with z_thresh={z_thresh}\n")
else:
    df_u = df_res[df_res["interaction"].isin(union_interactions)].copy()

    # --- ABS distance matrix (interaction x cell_type) ---
    mat_dist = (
        df_u.pivot_table(index="interaction", columns="cell_type", values="dist", aggfunc="max")
            .fillna(0.0)
    )
    mat_dist_c, _, _ = cluster_rows_and_cols(mat_dist, method="average", metric="correlation")
    save_heatmap(
        mat_dist_c,
        outfile=os.path.join(plots_dir, "heatmap_outlier_union_distance.png"),
        title="Outlier union heatmap (distance from trend line)",
        cbar_label="Distance to trend line"
    )
    mat_dist_c.to_csv(os.path.join(plots_dir, "heatmap_outlier_union_distance.csv"))

    # --- Binary matrix (interaction x cell_type): 0/1 ---
    df_u["outlier01"] = df_u["outlier"].astype("int8")
    mat_bin = (
        df_u.pivot_table(index="interaction", columns="cell_type", values="outlier01", aggfunc="max")
            .fillna(0).astype("int8")
    )
    mat_bin_c, _, _ = cluster_rows_and_cols(mat_bin, method="average", metric="jaccard")
    save_heatmap(
        mat_bin_c,
        outfile=os.path.join(plots_dir, "heatmap_outlier_union_binary.png"),
        title="Outlier union heatmap (binary)",
        cbar_label="Outlier (0/1)"
    )
    mat_bin_c.to_csv(os.path.join(plots_dir, "heatmap_outlier_union_binary.csv"))

    # --- SIGNED distance matrix (interaction x cell_type): dist_signed (neg/pos) ---
    mat_sdist = (
        df_u.pivot_table(index="interaction", columns="cell_type", values="dist_signed", aggfunc=max_abs_agg)
            .fillna(0.0)
    )
    mat_sdist_c, _, _ = cluster_rows_and_cols(mat_sdist, method="average", metric="correlation")
    save_heatmap(
        mat_sdist_c,
        outfile=os.path.join(plots_dir, "heatmap_outlier_union_signed_distance.png"),
        title="Outlier union heatmap (signed distance; direction)",
        cbar_label="Signed distance (neg=extracellular, pos=cellular)",
        cmap="coolwarm",
        center0=True
    )
    mat_sdist_c.to_csv(os.path.join(plots_dir, "heatmap_outlier_union_signed_distance.csv"))

    # --- SIGNED BINARY matrix (interaction x cell_type): -1/0/+1 (only for outliers) ---
    # 0 = not an outlier in this cell_type
    # +1 = outlier above line (residual > 0; cellular-leaning)
    # -1 = outlier below line (residual < 0; extracellular-leaning)
    df_u["signed_outlier"] = 0
    m = df_u["outlier"].to_numpy()
    df_u.loc[m, "signed_outlier"] = np.sign(df_u.loc[m, "residual"]).astype("int8")

    mat_sbin = (
        df_u.pivot_table(index="interaction", columns="cell_type", values="signed_outlier", aggfunc=max_abs_agg)
            .fillna(0)
            .astype("int8")
    )
    mat_sbin_c, _, _ = cluster_rows_and_cols(mat_sbin, method="average", metric="correlation")
    save_heatmap(
        mat_sbin_c,
        outfile=os.path.join(plots_dir, "heatmap_outlier_union_signed_binary.png"),
        title="Outlier union heatmap (signed binary; direction)",
        cbar_label="Signed outlier (neg=extracellular, pos=cellular; 0=not outlier)",
        cmap="coolwarm",
        center0=True
    )
    mat_sbin_c.to_csv(os.path.join(plots_dir, "heatmap_outlier_union_signed_binary.csv"))

    # --- Optional: save union outlier rows (long form) ---
    (
        df_u.sort_values(["cell_type", "dist"], ascending=[True, False])
           .to_csv(os.path.join(plots_dir, "outlier_union_table.csv"), index=False)
    )

print("Done.")
print(f"- Rank–rank plots: {rank_plot_dir}")
print(f"- Merged outlier table: {df_res_outfile}")
print(f"- Heatmaps saved under: {plots_dir}")
