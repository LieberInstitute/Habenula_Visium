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
import seaborn as sns

#   Read input files
in_dir = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana2'
)
plot_dir= here(
    'plots', '10_HD_bin_level', 'no_secondary', 'liana2'
)
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(os.path.join(in_dir, "table"), exist_ok=True)

in_files = [ os.path.join(in_dir, f) for f in os.listdir(in_dir) if re.compile(r'.*\.h5ad$').match(f) ]
in_files = [f for f in in_files if "extracellular" not in f]

#   Read in DataFrames of ligand-receptor stats for each donor and concatenate
lr_df_list = []
for f in in_files:
    donor_id = f.split('/')[-1].replace('.h5ad', '').replace('lrdata_', '')
    lr_df = sc.read(f).uns["global_interactions"]
    lr_df['donor_id'] = donor_id
    lr_df_list.append(lr_df)

lr_df = pd.concat(lr_df_list, axis=0, ignore_index=True)
lr_df = lr_df.rename(columns={"lr_mean": "mean", "ligand_complex": "ligand", "receptor_complex": "receptor"})

min_num_donors=10
#   Require a pair to be present in some minimum number of donors
lr_df = lr_df[
    lr_df.groupby(["source", "target", "ligand", "receptor"])["donor_id"]
    .transform("nunique") >= min_num_donors
].copy()

lr_df_sig = lr_df[
    (lr_df["pval"] < 0.05) &
    (lr_df["mean"] > 0)
].copy()

lr_df_sig.to_csv(os.path.join(in_dir, "table", "significant_interactions_across_donors.csv"), index=False)

#   Average stats across samples
lr_mean_df = (
    lr_df_sig
    .groupby(
        ["source", "target", "ligand", "receptor"],
        as_index=False,
        observed=True
    )
    .agg({"mean": "mean"})
)

lr_top_df = lr_mean_df.drop_duplicates(subset=['source', 'target', 'ligand', 'receptor'])
lr_top_df = lr_top_df.sort_values("mean", ascending=False)
lr_top_df.to_csv(os.path.join(in_dir, "table", "overall_mean_across_donors.csv"), index=False)
# ========================================

# -----------------------------
# figure1:top 20 by mean
plot_df = (
    lr_top_df
    .sort_values("mean", ascending=False)
    .head(20)
    .copy()
)

plot_df["interaction"] = (
    plot_df["source"].astype(str) + " → " + 
    plot_df["target"].astype(str) + "\n" +
    plot_df["ligand"].astype(str) + " - " + 
    plot_df["receptor"].astype(str)
)

top_palette1 = sns.color_palette("viridis", n_colors=len(plot_df))

fig, ax = plt.subplots(figsize=(14, 7))

sns.barplot(
    data=plot_df,
    x="interaction",
    y="mean",
    palette=top_palette1,
    ax=ax,
    edgecolor="black",
    linewidth=0.6,
)

ax.set_title("Top 20 Interactions by Mean", fontsize=18, weight="bold", pad=15)
ax.set_xlabel("")
ax.set_ylabel("Mean", fontsize=14)
ax.tick_params(axis="x", rotation=70, labelsize=10)
ax.tick_params(axis="y", labelsize=11)

for spine in ["top", "right"]:
    ax.spines[spine].set_visible(False)

plt.tight_layout()
plt.show()

plt.savefig(os.path.join(plot_dir, f"top20_mean_interactions.pdf"), dpi=300, bbox_inches="tight")
plt.close()