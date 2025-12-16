# Create the LIANA files for the multidonor analysis
import scanpy as sc
import decoupler as dc
import plotnine as p9
import liana as li
import numpy as np
import pandas as pd
import os
from pyhere import here
import matplotlib.pyplot as plt
import json
import numpy as np
from pyhere import here
import re
import session_info
from plotnine import *
from scipy import stats

#   Read input files
in_dir = here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'liana'
)
plot_dir= here(
    'processed-data', '10_HD_bin_level', 'new_samples2', 'liana', 'figure','habenula'
)
os.makedirs(plot_dir, exist_ok=True)

in_files = [ os.path.join(in_dir, f) for f in os.listdir(in_dir) if re.compile(r'.*\.h5ad$').match(f) ]
in_files = [f for f in in_files if "extracellular" not in f]

#   Read in DataFrames of ligand-receptor stats for each donor and concatenate
lr_df_list = []
for f in in_files:
    donor_id = f.split('/')[-1].replace('.h5ad', '').replace('lrdata_', '')
    lr_df = sc.read(f).var
    lr_df['donor_id'] = donor_id
    lr_df_list.append(lr_df)

lr_df = pd.concat(lr_df_list, axis=0, ignore_index=True)

min_num_donors=5
#   Require a pair to be present in some minimum number of donors
lr_df = lr_df[
    lr_df.groupby(['ligand', 'receptor'], observed=True)['ligand'].transform('count') >= min_num_donors
]

#   Average stats across samples
lr_mean_df = (
    lr_df
        .groupby(['ligand', 'receptor'], as_index=False)
        .agg({'mean': 'mean', 'morans': 'mean'})
)

lr_mean_df.to_csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples2/liana/table/overall_mean_morans_across_donors.csv", index=False)

lr_mean_df_morans = lr_mean_df.sort_values("morans", ascending=False)
lr_mean_df_mean = lr_mean_df.sort_values("mean", ascending=False)

lr_top_df = pd.concat(
    [
        lr_mean_df.sort_values("morans", ascending=False).head(),
        lr_mean_df.sort_values("mean", ascending=False).head()
    ]
)
lr_top_df = lr_top_df.drop_duplicates(subset=['ligand', 'receptor'])

top_pairs = list(lr_top_df['ligand'] + '^' + lr_top_df['receptor'])

#  first sample
for f in in_files:                      
    ad_lr = sc.read(f)                  
    sample_id = os.path.basename(f).replace(".h5ad", "").replace("lrdata_", "")

    assert all([x in ad_lr.var.index for x in top_pairs])

    #   Fix spatial coordinates format
    ad_lr.obsm['spatial'] = np.asarray(ad_lr.obsm['spatial'])

    #   Plot scores, permutation-based p-values, and local categories
    fig = sc.pl.spatial(
        ad_lr, color=top_pairs, cmap='viridis', ncols=2,
        spot_size=80, show=False, return_fig=True
    )
    fig.savefig(os.path.join(plot_dir, f"bivariate_scores_{sample_id}.pdf"), bbox_inches='tight')
    plt.close(fig)

    sc.pl.spatial(
        ad_lr, color=top_pairs, cmap='viridis_r', ncols=2,
        spot_size=80, layer='pvals'
    )
    plt.savefig(os.path.join(plot_dir, f"bivariate_p_vals_{sample_id}.pdf"))
    plt.close('all')

    sc.pl.spatial(
        ad_lr, color=top_pairs, cmap='coolwarm', ncols=2,
        spot_size=80, layer='cats'
    )
    plt.savefig(os.path.join(plot_dir, f"bivariate_cats_{sample_id}.pdf"))
    plt.close('all')

session_info.show()

# ========================================
# make figure based on the mean morans and mean mean
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

lr_mean_df["interaction"] = lr_mean_df["ligand"] + "^" + lr_mean_df["receptor"]
top_mean = lr_mean_df.sort_values("mean", ascending=False).head(20)
top_morans = lr_mean_df.sort_values("morans", ascending=False).head(20)

sns.set(style="whitegrid", context="talk")
top_palette1 = sns.color_palette("Blues_d", n_colors=20)
top_palette2 = sns.color_palette("Purples_d", n_colors=20)


fig, axes = plt.subplots(2, 1, figsize=(13, 11), sharex=False, gridspec_kw={"hspace": 1})
# -----------------------------
# figure1:top 20 by mean
sns.barplot(
    data=top_mean,
    x="interaction", y="mean",
    palette=top_palette1,
    ax=axes[0],
    edgecolor="black",
    linewidth=0.6,
)
axes[0].set_title("Top 20 Interactions by Mean", fontsize=18, weight="bold", pad=15)
axes[0].set_xlabel("")
axes[0].set_ylabel("Mean", fontsize=14)
axes[0].tick_params(axis="x", rotation=70, labelsize=10)
axes[0].tick_params(axis="y", labelsize=11)
for spine in ["top", "right"]:
    axes[0].spines[spine].set_visible(False)


# -----------------------------
# figure2: top 20 by morans
sns.barplot(
    data=top_morans,
    x="interaction", y="morans",
    palette=top_palette2,
    ax=axes[1],
    edgecolor="black",
    linewidth=0.6,
)
axes[1].set_title("Top 20 Interactions by Moran’s I", fontsize=18, weight="bold", pad=15)
axes[1].set_xlabel("Ligand–Receptor Pair", fontsize=14, labelpad=10)
axes[1].set_ylabel("Moran’s I", fontsize=14)
axes[1].tick_params(axis="x", rotation=70, labelsize=10)
axes[1].tick_params(axis="y", labelsize=11)
for spine in ["top", "right"]:
    axes[1].spines[spine].set_visible(False)

plt.tight_layout()
plt.savefig("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples2/liana/figure/habenula/Top20_LR_mean_morans.pdf", bbox_inches="tight")
plt.close()

print("✅ Figure saved as: Top20_LR_mean_morans.pdf")

# ---------------------------------
# scatter plot mean vs morans
plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=lr_mean_df,
    x="mean",
    y="morans",
    alpha=0.6,
    edgecolor=None,
)
plt.title("Ligand-Receptor Interactions: Mean vs Moran’s I", fontsize=18, weight="bold", pad=15)
plt.xlabel("Mean", fontsize=14)
plt.ylabel("Moran’s I", fontsize=14)
plt.grid(False)
plt.tight_layout()
plt.savefig("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples2/liana/figure/habenula/Mean_vs_Morans_scatter.pdf", bbox_inches="tight")
plt.close()
print("✅ Figure saved as: Mean_vs_Morans_scatter.pdf")

# ---------------------------------------
# scatter plot mean vs morans
slope, intercept, r_value, p_value, std_err = stats.linregress(lr_mean_df["mean"], lr_mean_df["morans"])
lr_mean_df["predicted"] = intercept + slope * lr_mean_df["mean"]
lr_mean_df["residual"] = lr_mean_df["morans"] - lr_mean_df["predicted"]

threshold = 2 * lr_mean_df["residual"].std()
outliers = lr_mean_df[np.abs(lr_mean_df["residual"]) > threshold]

plt.figure(figsize=(10, 8))
sns.scatterplot(
    data=lr_mean_df,
    x="mean",
    y="morans",
    alpha=0.6,
    edgecolor=None,
    label="Data"
)
sns.lineplot(
    x=lr_mean_df["mean"],
    y=lr_mean_df["predicted"],
    color="red",
    linewidth=2,
    label="Fitted line"
)

plt.scatter(outliers["mean"], outliers["morans"], color="orange", edgecolor="black", s=80, label="Outliers")

for _, row in outliers.iterrows():
    plt.text(
        row["mean"],
        row["morans"],
        row["interaction"],
        fontsize=9,
        color="black",
        ha="right",
        va="bottom"
    )

plt.title("Ligand–Receptor Interactions: Mean vs Moran’s I", fontsize=18, weight="bold", pad=15)
plt.xlabel("Mean", fontsize=14)
plt.ylabel("Moran’s I", fontsize=14)
plt.grid(False)
plt.legend()
plt.tight_layout()
plt.savefig("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/new_samples2/liana/figure/habenula/Mean_vs_Morans_scatter_outlier.pdf", bbox_inches="tight")
plt.close()
print("✅ Figure saved as: Mean_vs_Morans_scatter_outlier.pdf")