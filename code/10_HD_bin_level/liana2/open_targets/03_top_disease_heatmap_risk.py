import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pyhere import here

# Choose bandwidth
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--bandwidth", type=float, required=True)
parser.add_argument("--disease", type=str, required=True)
args = parser.parse_args()

bandwidth = args.bandwidth
disease_name = args.disease

print("bandwidth =", bandwidth)
print("disease =", disease_name)

plot_dir = here(
    'plots', '10_HD_bin_level', 'no_secondary', 'liana2',
    'open_targets'
)

out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana2',
    'open_targets'
)

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(out_path, exist_ok=True)

disease_pairs = pd.read_csv(
    os.path.join(out_path, f"{disease_name}_top_interactions.csv")
)

# --------------------------------- Fig 1 ---------------------------------

df = pd.read_csv(
    f"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/"
    f"processed-data/10_HD_bin_level/no_secondary/liana2/table/"
    f"significant_interactions_across_donors_{bandwidth}.csv"
)
df_before = df.copy()

df = df.merge(
    disease_pairs[["genesymbol_intercell_source", "genesymbol_intercell_target"]],
    left_on=["ligand", "receptor"],
    right_on=["genesymbol_intercell_source", "genesymbol_intercell_target"],
    how="inner"
).copy()

print("before:", df_before.shape)
print("after:", df.shape)

# -----------------------------
# Optional: filter for significant interactions first
# If you want to keep only significant ones, uncomment the next line
# df = df[df["pval"] < 0.05].copy()
# -----------------------------

# 1. First, count how many interactions each source-target pair has in each donor

count_per_donor = (
    df.groupby(["donor_id", "source", "target"])
      .size()
      .reset_index(name="n_interactions")
)

# 2. Then, average the number of interactions across donors for each source-target pair

avg_count = (
    count_per_donor.groupby(["source", "target"])["n_interactions"]
    .mean()
    .reset_index(name="avg_n_interactions")
)

# 3. transform to heatmap
heatmap_data = avg_count.pivot(
    index="target",     # y
    columns="source",   # x
    values="avg_n_interactions"
).fillna(0)

# 4. plot  heatmap

plt.figure(figsize=(11, 9))
sns.heatmap(
    heatmap_data,
    annot=True,
    fmt=".2f",
    cmap="YlOrRd",
    linewidths=0.5,
    linecolor="white"
)

plt.xlabel("Source")
plt.ylabel("Target")
plt.title("Average Number of Interactions per Donor")
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()

plt.savefig(    
    os.path.join(plot_dir,f"source_target_avgcount_heatmap_top_{disease_name}_{bandwidth}.png"), dpi=300, bbox_inches="tight")
plt.close()

# --------------------- Fig 2 ---------------------
heat_df = (
    df.groupby(["source", "target"], as_index=False)["mean"]
      .mean()
      .rename(columns={"mean": "avg_mean"})
)

heatmap_data = heat_df.pivot(
    index="target",      # y
    columns="source",    # x
    values="avg_mean"
).fillna(0)

plt.figure(figsize=(11, 9))
sns.heatmap(
    heatmap_data,
    annot=True,
    fmt=".2f",
    cmap="YlGnBu",
    linewidths=0.5,
    linecolor="white"
)

plt.xlabel("Source")
plt.ylabel("Target")
plt.title("Average LR Interaction Mean by Source-Target")
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()

plt.savefig(os.path.join(plot_dir,f"source_target_mean_heatmap_top_{disease_name}_{bandwidth}.png") , dpi=300, bbox_inches="tight")
plt.close()


# --------------------- Fig 3 ---------------------
heat_df = (
    df.groupby(["source", "target"], as_index=False)["mean"]
      .sum()
      .rename(columns={"mean": "sum_mean"})
)

heatmap_data = heat_df.pivot(
    index="target",      # y
    columns="source",    # x
    values="sum_mean"
).fillna(0)

plt.figure(figsize=(11, 9))
sns.heatmap(
    heatmap_data,
    annot=True,
    fmt=".2f",
    cmap="YlOrRd",
    linewidths=0.5,
    linecolor="white"
)

plt.xlabel("Source")
plt.ylabel("Target")
plt.title("Summed LR Interaction Mean by Source-Target")
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()

plt.savefig(os.path.join(plot_dir,f"source_target_sum_heatmap_top_{disease_name}_{bandwidth}.png"), dpi=300, bbox_inches="tight")
plt.close()


