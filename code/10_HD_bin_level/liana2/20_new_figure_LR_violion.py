# bandwidth vs. vilion plot
# make a figure

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
import re
from pyhere import here

# -----------------------------
# Paths
# -----------------------------
table_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana2/table"

plot_dir = here(
    "plots", "10_HD_bin_level", "no_secondary", "liana2"
)
os.makedirs(plot_dir, exist_ok=True)

# Whether to include source-target pairs that are missing in a bandwidth as 0
# Usually False is better, because fillna(0) creates artificial zeros.
INCLUDE_MISSING_ZEROS = False

# -----------------------------
# Helper function
# -----------------------------
def extract_bandwidth(filename):
    """
    Extract bandwidth from file name.
    Example:
    significant_interactions_across_donors_2500.0.csv -> 2500.0
    """
    match = re.search(r"_(\d+\.?\d*)\.csv$", filename)
    if match:
        return float(match.group(1))
    return None


# -----------------------------
# Part 1: Average number of interactions
# from significant_interactions_across_donors_{bandwidth}.csv
# -----------------------------
avg_count_list = []

sig_files = glob.glob(
    os.path.join(table_dir, "significant_interactions_across_donors_*.csv")
)

for file in sig_files:
    bandwidth = extract_bandwidth(os.path.basename(file))
    if bandwidth is None:
        continue

    df = pd.read_csv(file)

    count_per_donor = (
        df.groupby(["donor_id", "source", "target"])
          .size()
          .reset_index(name="n_interactions")
    )

    avg_count = (
        count_per_donor.groupby(["source", "target"])["n_interactions"]
        .mean()
        .reset_index(name="value")
    )

    if INCLUDE_MISSING_ZEROS:
        heatmap_data = avg_count.pivot(
            index="target",
            columns="source",
            values="value"
        ).fillna(0)

        values = heatmap_data.values.flatten()

        temp = pd.DataFrame({
            "bandwidth": bandwidth,
            "value": values
        })

    else:
        temp = avg_count[["value"]].copy()
        temp["bandwidth"] = bandwidth

    temp["metric"] = "Average number of interactions"
    avg_count_list.append(temp)


avg_count_violin = pd.concat(avg_count_list, ignore_index=True)


# -----------------------------
# Part 2: Average LR mean and summed LR mean
# from overall_mean_across_donors_{bandwidth}.csv
# -----------------------------
avg_mean_list = []
sum_mean_list = []

mean_files = glob.glob(
    os.path.join(table_dir, "overall_mean_across_donors_*.csv")
)

for file in mean_files:
    bandwidth = extract_bandwidth(os.path.basename(file))
    if bandwidth is None:
        continue

    df = pd.read_csv(file)

    # Average mean per source-target pair
    avg_mean = (
        df.groupby(["source", "target"], as_index=False)["mean"]
          .mean()
          .rename(columns={"mean": "value"})
    )

    if INCLUDE_MISSING_ZEROS:
        heatmap_data = avg_mean.pivot(
            index="target",
            columns="source",
            values="value"
        ).fillna(0)

        values = heatmap_data.values.flatten()

        temp_avg = pd.DataFrame({
            "bandwidth": bandwidth,
            "value": values
        })

    else:
        temp_avg = avg_mean[["value"]].copy()
        temp_avg["bandwidth"] = bandwidth

    temp_avg["metric"] = "Average LR interaction mean"
    avg_mean_list.append(temp_avg)

    # Summed mean per source-target pair
    sum_mean = (
        df.groupby(["source", "target"], as_index=False)["mean"]
          .sum()
          .rename(columns={"mean": "value"})
    )

    if INCLUDE_MISSING_ZEROS:
        heatmap_data = sum_mean.pivot(
            index="target",
            columns="source",
            values="value"
        ).fillna(0)

        values = heatmap_data.values.flatten()

        temp_sum = pd.DataFrame({
            "bandwidth": bandwidth,
            "value": values
        })

    else:
        temp_sum = sum_mean[["value"]].copy()
        temp_sum["bandwidth"] = bandwidth

    temp_sum["metric"] = "Summed LR interaction mean"
    sum_mean_list.append(temp_sum)


avg_mean_violin = pd.concat(avg_mean_list, ignore_index=True)
sum_mean_violin = pd.concat(sum_mean_list, ignore_index=True)


# -----------------------------
# Plot function
# -----------------------------
def plot_violin(data, y_label, title, output_name):
    data = data.copy()
    data["bandwidth"] = data["bandwidth"].astype(str)

    # Sort bandwidth numerically
    order = (
        data[["bandwidth"]]
        .drop_duplicates()
        .assign(bw_num=lambda x: x["bandwidth"].astype(float))
        .sort_values("bw_num")["bandwidth"]
        .tolist()
    )

    plt.figure(figsize=(10, 6))

    sns.violinplot(
        data=data,
        x="bandwidth",
        y="value",
        order=order,
        inner=None,
        cut=0
    )

    sns.stripplot(
        data=data,
        x="bandwidth",
        y="value",
        order=order,
        jitter=True,
        size=3,
        alpha=0.5
    )

    plt.xlabel("Bandwidth")
    plt.ylabel(y_label)
    plt.title(title)
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    plt.savefig(
        os.path.join(plot_dir, output_name),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


# -----------------------------
# Make plots
# -----------------------------
plot_violin(
    avg_count_violin,
    y_label="Average number of interactions",
    title="Distribution of Heatmap Values Across Bandwidths",
    output_name="violin_avg_n_interactions_by_bandwidth.png"
)

plot_violin(
    avg_mean_violin,
    y_label="Average LR interaction mean",
    title="Distribution of Average LR Mean Across Bandwidths",
    output_name="violin_avg_lr_mean_by_bandwidth.png"
)

plot_violin(
    sum_mean_violin,
    y_label="Summed LR interaction mean",
    title="Distribution of Summed LR Mean Across Bandwidths",
    output_name="violin_sum_lr_mean_by_bandwidth.png"
)