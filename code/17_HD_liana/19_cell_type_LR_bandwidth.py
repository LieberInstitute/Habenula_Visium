#!/usr/bin/env python3

import argparse
import os

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pyhere import here


###############################################################################
# User settings
###############################################################################

CELL_TYPE_MAP_FILE = (
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/"
    "Habenula_Visium/raw-data/sample_info/"
    "hd_cell_type_map.csv"
)

PROCESSED_DATA_DIR = (
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/"
    "Habenula_Visium/processed-data/10_HD_bin_level/"
    "no_secondary/liana2/table"
)

PLOT_DIR = here(
    "plots",
    "10_HD_bin_level",
    "no_secondary",
    "liana2"
)


###############################################################################
# Command-line arguments
###############################################################################

parser = argparse.ArgumentParser(
    description=(
        "Plot LIANA source-target heatmaps using corrected cell-type names."
    )
)

parser.add_argument(
    "--bandwidth",
    type=float,
    required=True,
    help="Bandwidth value used in the LIANA result filenames."
)

args = parser.parse_args()
bandwidth = args.bandwidth

print(f"bandwidth = {bandwidth}")

os.makedirs(PLOT_DIR, exist_ok=True)


###############################################################################
# Helper functions
###############################################################################

def read_cell_type_mapping(mapping_file):
    """
    Read and validate the old-to-new cell-type name mapping.
    """

    if not os.path.exists(mapping_file):
        raise FileNotFoundError(
            f"Cell-type mapping file does not exist: {mapping_file}"
        )

    if os.path.getsize(mapping_file) == 0:
        raise ValueError(
            f"Cell-type mapping file is empty: {mapping_file}"
        )

    mapping_df = pd.read_csv(
        mapping_file,
        dtype=str
    )

    required_columns = {
        "old_cell_type",
        "new_cell_type",
        "color"
    }

    missing_columns = required_columns.difference(mapping_df.columns)

    if missing_columns:
        raise ValueError(
            "Missing required columns in cell-type mapping file: "
            + ", ".join(sorted(missing_columns))
        )

    for column in required_columns:
        mapping_df[column] = mapping_df[column].str.strip()

    mapping_df = mapping_df.dropna(
        subset=[
            "old_cell_type",
            "new_cell_type"
        ]
    )

    mapping_df = mapping_df[
        (mapping_df["old_cell_type"] != "") &
        (mapping_df["new_cell_type"] != "")
    ].copy()

    if mapping_df.empty:
        raise ValueError(
            "No valid rows remained in the cell-type mapping file."
        )

    duplicated_old_names = mapping_df.loc[
        mapping_df["old_cell_type"].duplicated(keep=False),
        [
            "old_cell_type",
            "new_cell_type",
            "color"
        ]
    ]

    if not duplicated_old_names.empty:
        print("\nDuplicated old_cell_type entries:")
        print(
            duplicated_old_names
            .sort_values("old_cell_type")
            .to_string(index=False)
        )

        raise ValueError(
            "Each old_cell_type must appear only once in the mapping file."
        )

    name_map = dict(
        zip(
            mapping_df["old_cell_type"],
            mapping_df["new_cell_type"]
        )
    )

    color_map = dict(
        zip(
            mapping_df["new_cell_type"],
            mapping_df["color"]
        )
    )

    return mapping_df, name_map, color_map


def rename_source_target(df, name_map, data_name):
    """
    Replace source and target names using the provided mapping.

    Unmatched names are retained and reported.
    """

    required_columns = {
        "source",
        "target"
    }

    missing_columns = required_columns.difference(df.columns)

    if missing_columns:
        raise ValueError(
            f"{data_name} is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    df = df.copy()

    df["source"] = df["source"].astype(str).str.strip()
    df["target"] = df["target"].astype(str).str.strip()

    # Preserve the original names for auditing.
    df["source_original"] = df["source"]
    df["target_original"] = df["target"]

    known_old_names = set(name_map)

    unmatched_sources = sorted(
        set(df["source"]) - known_old_names
    )

    unmatched_targets = sorted(
        set(df["target"]) - known_old_names
    )

    # Replace matched names and retain unmatched names.
    df["source"] = df["source"].map(name_map).fillna(
        df["source_original"]
    )

    df["target"] = df["target"].map(name_map).fillna(
        df["target_original"]
    )

    print(f"\nName mapping summary for {data_name}:")
    print(
        "  Unique source names: "
        f"{df['source_original'].nunique()}"
    )
    print(
        "  Unique target names: "
        f"{df['target_original'].nunique()}"
    )
    print(
        "  Mapped source names: "
        f"{df.loc[df['source_original'].isin(known_old_names), 'source_original'].nunique()}"
    )
    print(
        "  Mapped target names: "
        f"{df.loc[df['target_original'].isin(known_old_names), 'target_original'].nunique()}"
    )

    if unmatched_sources:
        print(
            "\nWARNING: These source names were not found in the mapping "
            "file and will retain their original names:"
        )
        for name in unmatched_sources:
            print(f"  - {name}")

    if unmatched_targets:
        print(
            "\nWARNING: These target names were not found in the mapping "
            "file and will retain their original names:"
        )
        for name in unmatched_targets:
            print(f"  - {name}")

    return df


def order_heatmap_axes(
    heatmap_data,
    cell_type_order
):
    """
    Order heatmap columns and rows according to the mapping-file order.

    Any unmatched names are placed after the mapped cell types.
    """

    mapped_columns = [
        cell_type
        for cell_type in cell_type_order
        if cell_type in heatmap_data.columns
    ]

    extra_columns = sorted(
        set(heatmap_data.columns) - set(mapped_columns)
    )

    mapped_rows = [
        cell_type
        for cell_type in cell_type_order
        if cell_type in heatmap_data.index
    ]

    extra_rows = sorted(
        set(heatmap_data.index) - set(mapped_rows)
    )

    return heatmap_data.reindex(
        index=mapped_rows + extra_rows,
        columns=mapped_columns + extra_columns,
        fill_value=0
    )


def save_heatmap(
    heatmap_data,
    output_file,
    title,
    cmap,
    fmt=".1f"
):
    """
    Plot and save a source-target heatmap.
    """

    n_sources = heatmap_data.shape[1]
    n_targets = heatmap_data.shape[0]

    figure_width = max(
        11,
        0.75 * n_sources + 3
    )

    figure_height = max(
        9,
        0.55 * n_targets + 3
    )

    plt.figure(
        figsize=(
            figure_width,
            figure_height
        )
    )

    ax = sns.heatmap(
        heatmap_data,
        annot=True,
        fmt=fmt,
        cmap=cmap,
        linewidths=0.5,
        linecolor="white",
        cbar=True,
        annot_kws={
            "fontsize": 12
        }
    )

    # Axis labels
    ax.set_xlabel(
        "Source",
        fontsize=16
    )

    ax.set_ylabel(
        "Target",
        fontsize=16
    )

    # Title
    ax.set_title(
        title,
        fontsize=18,
        pad=15
    )

    # X-axis tick labels
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        fontsize=13
    )

    # Y-axis tick labels
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=13
    )

    # Colorbar tick labels
    cbar = ax.collections[0].colorbar
    cbar.ax.tick_params(
        labelsize=13
    )

    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Saved: {output_file}")


###############################################################################
# Read cell-type mapping
###############################################################################

cell_type_map_df, cell_type_name_map, cell_type_color_map = (
    read_cell_type_mapping(
        CELL_TYPE_MAP_FILE
    )
)

# Use the mapping-file order for heatmap rows and columns.
cell_type_order = (
    cell_type_map_df["new_cell_type"]
    .drop_duplicates()
    .tolist()
)

print("\nCell-type mapping:")
print(
    cell_type_map_df[
        [
            "old_cell_type",
            "new_cell_type",
            "color"
        ]
    ].to_string(index=False)
)


###############################################################################
# Figure 1: Average number of interactions per donor
###############################################################################

significant_interactions_file = os.path.join(
    PROCESSED_DATA_DIR,
    f"significant_interactions_across_donors_{bandwidth}.csv"
)

if not os.path.exists(significant_interactions_file):
    raise FileNotFoundError(
        "Significant-interactions file does not exist: "
        f"{significant_interactions_file}"
    )

significant_df = pd.read_csv(
    significant_interactions_file
)

significant_df = rename_source_target(
    significant_df,
    cell_type_name_map,
    data_name=os.path.basename(significant_interactions_file)
)

# Count interactions for each donor and source-target pair.
count_per_donor = (
    significant_df
    .groupby(
        [
            "donor_id",
            "source",
            "target"
        ],
        as_index=False
    )
    .size()
    .rename(
        columns={
            "size": "n_interactions"
        }
    )
)

# Average interaction count across donors.
avg_count = (
    count_per_donor
    .groupby(
        [
            "source",
            "target"
        ],
        as_index=False
    )["n_interactions"]
    .mean()
    .rename(
        columns={
            "n_interactions": "avg_n_interactions"
        }
    )
)

heatmap_avg_count = (
    avg_count
    .pivot(
        index="target",
        columns="source",
        values="avg_n_interactions"
    )
    .fillna(0)
)

heatmap_avg_count = order_heatmap_axes(
    heatmap_avg_count,
    cell_type_order
)

avg_count_data_file = os.path.join(
    PLOT_DIR,
    f"source_target_avgcount_heatmap_data_{bandwidth}.csv"
)

heatmap_avg_count.to_csv(
    avg_count_data_file,
    index=True
)

print(f"Saved: {avg_count_data_file}")

avg_count_plot_file = os.path.join(
    PLOT_DIR,
    f"source_target_avgcount_heatmap_{bandwidth}.pdf"
)

save_heatmap(
    heatmap_data=heatmap_avg_count,
    output_file=avg_count_plot_file,
    title="Average Number of Interactions per Donor",
    cmap="YlOrRd",
    fmt=".1f"
)


###############################################################################
# Read data for Figures 2 and 3
###############################################################################

overall_mean_file = os.path.join(
    PROCESSED_DATA_DIR,
    f"overall_mean_across_donors_{bandwidth}.csv"
)

if not os.path.exists(overall_mean_file):
    raise FileNotFoundError(
        "Overall-mean file does not exist: "
        f"{overall_mean_file}"
    )

mean_df = pd.read_csv(
    overall_mean_file
)

required_mean_columns = {
    "source",
    "target",
    "mean"
}

missing_mean_columns = required_mean_columns.difference(
    mean_df.columns
)

if missing_mean_columns:
    raise ValueError(
        "Overall-mean file is missing required columns: "
        + ", ".join(sorted(missing_mean_columns))
    )

mean_df["mean"] = pd.to_numeric(
    mean_df["mean"],
    errors="coerce"
)

n_before = len(mean_df)

mean_df = mean_df.dropna(
    subset=["mean"]
).copy()

n_removed = n_before - len(mean_df)

if n_removed > 0:
    print(
        f"\nRemoved {n_removed} rows with missing or invalid mean values."
    )

mean_df = rename_source_target(
    mean_df,
    cell_type_name_map,
    data_name=os.path.basename(overall_mean_file)
)


###############################################################################
# Figure 2: Average LR interaction mean
###############################################################################

heat_mean_df = (
    mean_df
    .groupby(
        [
            "source",
            "target"
        ],
        as_index=False
    )["mean"]
    .mean()
    .rename(
        columns={
            "mean": "avg_mean"
        }
    )
)

heatmap_avg_mean = (
    heat_mean_df
    .pivot(
        index="target",
        columns="source",
        values="avg_mean"
    )
    .fillna(0)
)

heatmap_avg_mean = order_heatmap_axes(
    heatmap_avg_mean,
    cell_type_order
)

avg_mean_data_file = os.path.join(
    PLOT_DIR,
    f"source_target_mean_heatmap_data_{bandwidth}.csv"
)

heatmap_avg_mean.to_csv(
    avg_mean_data_file,
    index=True
)

print(f"Saved: {avg_mean_data_file}")

avg_mean_plot_file = os.path.join(
    PLOT_DIR,
    f"source_target_mean_heatmap_{bandwidth}.pdf"
)

save_heatmap(
    heatmap_data=heatmap_avg_mean,
    output_file=avg_mean_plot_file,
    title="Average LR Interaction Mean by Source-Target",
    cmap="YlGnBu",
    fmt=".1f"
)


###############################################################################
# Figure 3: Summed LR interaction mean
###############################################################################

heat_sum_df = (
    mean_df
    .groupby(
        [
            "source",
            "target"
        ],
        as_index=False
    )["mean"]
    .sum()
    .rename(
        columns={
            "mean": "sum_mean"
        }
    )
)

heatmap_sum_mean = (
    heat_sum_df
    .pivot(
        index="target",
        columns="source",
        values="sum_mean"
    )
    .fillna(0)
)

heatmap_sum_mean = order_heatmap_axes(
    heatmap_sum_mean,
    cell_type_order
)

sum_mean_data_file = os.path.join(
    PLOT_DIR,
    f"source_target_sum_mean_heatmap_data_{bandwidth}.csv"
)

heatmap_sum_mean.to_csv(
    sum_mean_data_file,
    index=True
)

print(f"Saved: {sum_mean_data_file}")

sum_mean_plot_file = os.path.join(
    PLOT_DIR,
    f"source_target_sum_heatmap_{bandwidth}.pdf"
)

save_heatmap(
    heatmap_data=heatmap_sum_mean,
    output_file=sum_mean_plot_file,
    title="Summed LR Interaction Mean by Source-Target",
    cmap="YlOrRd",
    fmt=".1f"
)


###############################################################################
# Save mapping audit
###############################################################################

mapping_audit_file = os.path.join(
    PLOT_DIR,
    f"cell_type_name_mapping_{bandwidth}.csv"
)

cell_type_map_df.to_csv(
    mapping_audit_file,
    index=False
)

print(f"Saved: {mapping_audit_file}")

print("\nAll heatmaps completed successfully.")