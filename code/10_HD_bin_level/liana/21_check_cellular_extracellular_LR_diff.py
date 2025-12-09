# Read in multi-donor data and plot top pairs
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# -----------------------
# 1) Load and merge (cellular vs extracellular)
# -----------------------
plots = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/figure"

cell_path = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/overall_mean_morans_across_donors.csv"
extra_path = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples/liana/table/overall_mean_morans_across_donors_extracellular.csv"

df_cell = pd.read_csv(cell_path)
df_extra = pd.read_csv(extra_path)

# Create a unique LR identifier
df_cell["LR"] = df_cell["ligand"].astype(str) + "_" + df_cell["receptor"].astype(str)
df_extra["LR"] = df_extra["ligand"].astype(str) + "_" + df_extra["receptor"].astype(str)

# Merge so each LR has both cellular and extracellular values
df = df_cell[["LR", "mean", "morans"]].merge(
    df_extra[["LR", "mean", "morans"]],
    on="LR",
    how="inner",
    suffixes=("_cellular", "_extracellular")
).dropna()

# -----------------------
# 2) Regression-based outliers using internal studentized residuals
#    Outlier rule: |studentized residual| >= threshold
# -----------------------
def regression_studentized_outliers(x, y, threshold=3.0):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n < 5:
        raise ValueError("Sample size is too small for stable regression/outlier detection.")

    # Design matrix X = [1, x]
    X = np.column_stack([np.ones(n), x])
    XtX_inv = np.linalg.inv(X.T @ X)
    beta = XtX_inv @ (X.T @ y)  # [intercept, slope]
    yhat = X @ beta
    resid = y - yhat

    # Leverage h_ii = diag(H), where H = X (X'X)^-1 X'
    H = X @ XtX_inv @ X.T
    h = np.clip(np.diag(H), 1e-12, 1 - 1e-12)

    # Mean squared error
    p = 2
    mse = (resid @ resid) / max(n - p, 1)

    # Internal studentized residuals
    stud = resid / np.sqrt(mse * (1 - h))

    out = np.abs(stud) >= threshold
    return {
        "intercept": beta[0],
        "slope": beta[1],
        "yhat": yhat,
        "resid": resid,
        "stud_resid": stud,
        "outlier": out
    }

# -----------------------
# 3) Make two plots: Pearson (raw values) & Spearman (rank-rank)
# -----------------------
def plot_pearson_and_spearman(df, metric="morans", outdir="plots", threshold=3.0, label_outliers=True):
    os.makedirs(outdir, exist_ok=True)

    xcol = f"{metric}_cellular"
    ycol = f"{metric}_extracellular"

    d = df[["LR", xcol, ycol]].dropna().copy()
    d = d[np.isfinite(d[xcol]) & np.isfinite(d[ycol])].copy()

    # ===== Pearson plot (raw values) =====
    pearson_r, pearson_p = stats.pearsonr(d[xcol].values, d[ycol].values)
    pear = regression_studentized_outliers(d[xcol].values, d[ycol].values, threshold=threshold)
    d["outlier_pearson"] = pear["outlier"]

    plt.figure(figsize=(7, 6))
    sns.scatterplot(
        data=d,
        x=xcol, y=ycol,
        hue="outlier_pearson",
        palette={False: "gray", True: "red"},
        s=40
    )

    # y = x reference line
    lims = [
        min(d[xcol].min(), d[ycol].min()),
        max(d[xcol].max(), d[ycol].max())
    ]
    plt.plot(lims, lims, "k--", linewidth=1)

    # Regression line
    xs = np.linspace(lims[0], lims[1], 200)
    ys = pear["intercept"] + pear["slope"] * xs
    plt.plot(xs, ys, linewidth=1)

    # Label outliers
    if label_outliers:
        for _, row in d[d["outlier_pearson"]].iterrows():
            plt.text(row[xcol], row[ycol], row["LR"], fontsize=7, color="red")

    plt.xlim(lims); plt.ylim(lims)
    plt.xlabel(f"{metric} (Cellular)")
    plt.ylabel(f"{metric} (Extracellular)")
    plt.title(f"Pearson (raw): r={pearson_r:.3f}, p={pearson_p:.2e}  |  outlier: |studentized resid|≥{threshold}")
    plt.tight_layout()

    plt.savefig(f"{outdir}/{metric}_pearson_raw.png", dpi=300)
    plt.savefig(f"{outdir}/{metric}_pearson_raw.pdf")
    plt.show()

    # ===== Spearman plot (rank-rank) =====
    # Spearman correlation is Pearson correlation on ranks
    d["rx"] = d[xcol].rank(method="average")
    d["ry"] = d[ycol].rank(method="average")

    spearman_rho, spearman_p = stats.spearmanr(d[xcol].values, d[ycol].values)
    spear = regression_studentized_outliers(d["rx"].values, d["ry"].values, threshold=threshold)
    d["outlier_spearman"] = spear["outlier"]

    plt.figure(figsize=(7, 6))
    sns.scatterplot(
        data=d,
        x="rx", y="ry",
        hue="outlier_spearman",
        palette={False: "gray", True: "red"},
        s=40
    )

    # y = x reference line in rank space
    r_lims = [
        min(d["rx"].min(), d["ry"].min()),
        max(d["rx"].max(), d["ry"].max())
    ]
    plt.plot(r_lims, r_lims, "k--", linewidth=1)

    # Regression line in rank space
    rxs = np.linspace(r_lims[0], r_lims[1], 200)
    rys = spear["intercept"] + spear["slope"] * rxs
    plt.plot(rxs, rys, linewidth=1)

    # Label outliers (in rank space)
    if label_outliers:
        for _, row in d[d["outlier_spearman"]].iterrows():
            plt.text(row["rx"], row["ry"], row["LR"], fontsize=7, color="red")

    plt.xlim(r_lims); plt.ylim(r_lims)
    plt.xlabel(f"Rank of {metric} (Cellular)")
    plt.ylabel(f"Rank of {metric} (Extracellular)")
    plt.title(f"Spearman (rank-rank): ρ={spearman_rho:.3f}, p={spearman_p:.2e}  |  outlier: |studentized resid|≥{threshold}")
    plt.tight_layout()

    plt.savefig(f"{outdir}/{metric}_spearman_rankrank.png", dpi=300)
    plt.savefig(f"{outdir}/{metric}_spearman_rankrank.pdf")
    plt.show()

# -----------------------
# 4) Run: generate plots for both mean and morans
# -----------------------
plot_pearson_and_spearman(df, metric="mean",   outdir=plots, threshold=2.0)
plot_pearson_and_spearman(df, metric="morans", outdir=plots, threshold=2.0)
