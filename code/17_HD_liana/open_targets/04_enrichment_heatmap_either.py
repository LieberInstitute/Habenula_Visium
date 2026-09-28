# Using the disease-risk data from opentargets, ask if the ligands/receptors
# are enriched for disease-associated genes within pairs of cell types.
# Heatmap axes = source and target cell types.
# Fill = enrichment statistic.
# Text = FDR significance stars only.

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import argparse

from pyhere import here
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


# -----------------------------
# args
# -----------------------------
parser = argparse.ArgumentParser()

parser.add_argument(
    "--bandwidth",
    type=float,
    required=True
)

parser.add_argument(
    "--disease",
    type=str,
    required=True
)

# Default: require LR interaction to be significant in all 8 donors
parser.add_argument(
    "--min_donors",
    type=int,
    default=8
)

args = parser.parse_args()

bandwidth = args.bandwidth
disease_name = args.disease
MIN_DONORS = args.min_donors


# -----------------------------
# Latest HD cell-type naming map
# -----------------------------
CELL_TYPE_MAP_FILE = (
    "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/"
    "Habenula_Visium/raw-data/sample_info/"
    "hd_cell_type_map.csv"
)


print("bandwidth =", bandwidth)
print("disease =", disease_name)
print("MIN_DONORS =", MIN_DONORS)


# -----------------------------
# paths
# -----------------------------
plot_dir = here(
    "plots",
    "10_HD_bin_level",
    "no_secondary",
    "liana2",
    "open_targets"
)

out_path = here(
    "processed-data",
    "10_HD_bin_level",
    "no_secondary",
    "liana2",
    "open_targets"
)

os.makedirs(
    plot_dir,
    exist_ok=True
)

os.makedirs(
    out_path,
    exist_ok=True
)


# -----------------------------
# load disease-risk ligand/receptor genes from OpenTargets
# -----------------------------
disease_pairs = pd.read_csv(
    os.path.join(
        out_path,
        f"{disease_name}_LorR_interactions.csv"
    )
)

# expected columns:
# genesymbol_intercell_source, genesymbol_intercell_target

disease_ligands_raw = set(
    disease_pairs[
        "genesymbol_intercell_source"
    ]
    .dropna()
    .unique()
)

disease_receptors_raw = set(
    disease_pairs[
        "genesymbol_intercell_target"
    ]
    .dropna()
    .unique()
)

print(
    "n raw disease ligands:",
    len(disease_ligands_raw)
)

print(
    "n raw disease receptors:",
    len(disease_receptors_raw)
)


# -----------------------------
# load LIANA interactions
# -----------------------------
df = pd.read_csv(
    f"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/"
    f"processed-data/10_HD_bin_level/no_secondary/liana2/table/"
    f"significant_interactions_across_donors_{bandwidth}.csv"
)


# -----------------------------
# map LIANA source/target to latest HD cell-type names
# -----------------------------
cell_type_map = pd.read_csv(
    CELL_TYPE_MAP_FILE
)

required_map_cols = {
    "old_cell_type",
    "new_cell_type"
}

missing_map_cols = (
    required_map_cols -
    set(cell_type_map.columns)
)

if missing_map_cols:
    raise ValueError(
        "CELL_TYPE_MAP_FILE is missing required columns: "
        f"{sorted(missing_map_cols)}"
    )


cell_type_map = (
    cell_type_map[
        [
            "old_cell_type",
            "new_cell_type"
        ]
    ]
    .dropna()
    .copy()
)

cell_type_map["old_cell_type"] = (
    cell_type_map["old_cell_type"]
    .astype(str)
    .str.strip()
)

cell_type_map["new_cell_type"] = (
    cell_type_map["new_cell_type"]
    .astype(str)
    .str.strip()
)


# Each old label should map to exactly one new label
n_map_per_old = (
    cell_type_map
    .groupby("old_cell_type")["new_cell_type"]
    .nunique()
)

ambiguous_old = n_map_per_old[
    n_map_per_old > 1
]

if not ambiguous_old.empty:
    raise ValueError(
        "Some old_cell_type values map to multiple "
        "new_cell_type values:\n"
        + ambiguous_old.to_string()
    )


cell_type_map = (
    cell_type_map
    .drop_duplicates("old_cell_type")
)

cell_type_dict = dict(
    zip(
        cell_type_map["old_cell_type"],
        cell_type_map["new_cell_type"]
    )
)


# Save originals for auditing;
# keep unmatched labels unchanged
df["source_original"] = df["source"]
df["target_original"] = df["target"]

df["source"] = (
    df["source"]
    .map(cell_type_dict)
    .fillna(df["source_original"])
)

df["target"] = (
    df["target"]
    .map(cell_type_dict)
    .fillna(df["target_original"])
)


print(
    "CELL_TYPE_MAP_FILE =",
    CELL_TYPE_MAP_FILE
)

print(
    "mapped source labels:",
    (
        df["source"] !=
        df["source_original"]
    ).sum(),
    "/",
    len(df)
)

print(
    "mapped target labels:",
    (
        df["target"] !=
        df["target_original"]
    ).sum(),
    "/",
    len(df)
)


# Use mapping-file order for heatmap axes,
# followed by any unmapped labels
mapped_cell_type_order = list(
    dict.fromkeys(
        cell_type_map[
            "new_cell_type"
        ].tolist()
    )
)

observed_cell_types = list(
    dict.fromkeys(
        pd.concat(
            [
                df["source"],
                df["target"]
            ],
            ignore_index=True
        ).tolist()
    )
)

CELL_TYPE_ORDER = (
    [
        x
        for x in mapped_cell_type_order
        if x in observed_cell_types
    ]
    +
    [
        x
        for x in observed_cell_types
        if x not in mapped_cell_type_order
    ]
)


# keep only needed columns
df = (
    df[
        [
            "donor_id",
            "source",
            "target",
            "ligand",
            "receptor"
        ]
    ]
    .dropna()
    .copy()
)

print(
    "input interaction rows:",
    df.shape[0]
)

print(
    "n donors:",
    df["donor_id"].nunique()
)


# -----------------------------
# filter LR interactions by donor support
# -----------------------------
donor_support = (
    df
    .groupby(
        [
            "source",
            "target",
            "ligand",
            "receptor"
        ]
    )["donor_id"]
    .nunique()
    .reset_index(
        name="n_donors"
    )
)

print(
    "unique source-target-LR before donor filter:",
    donor_support.shape[0]
)


df = df.merge(
    donor_support,
    on=[
        "source",
        "target",
        "ligand",
        "receptor"
    ],
    how="left"
)

df = df[
    df["n_donors"] >= MIN_DONORS
].copy()


print(
    f"interaction rows after donor filter "
    f"n_donors >= {MIN_DONORS}:",
    df.shape[0]
)

print(
    f"unique source-target-LR after donor filter "
    f"n_donors >= {MIN_DONORS}:",
    df[
        [
            "source",
            "target",
            "ligand",
            "receptor"
        ]
    ]
    .drop_duplicates()
    .shape[0]
)


# -----------------------------
# define LR-pair universe after donor filter
# -----------------------------
df_lr = (
    df[
        [
            "source",
            "target",
            "ligand",
            "receptor"
        ]
    ]
    .drop_duplicates()
    .copy()
)


# observable ligand/receptor universe
all_ligands = set(
    df_lr["ligand"].unique()
)

all_receptors = set(
    df_lr["receptor"].unique()
)


# restrict disease genes to genes
# observable in LIANA universe
disease_ligands = (
    disease_ligands_raw &
    all_ligands
)

disease_receptors = (
    disease_receptors_raw &
    all_receptors
)


print(
    "n disease ligands in LIANA universe:",
    len(disease_ligands)
)

print(
    "n disease receptors in LIANA universe:",
    len(disease_receptors)
)

print(
    "n total ligand universe:",
    len(all_ligands)
)

print(
    "n total receptor universe:",
    len(all_receptors)
)


# unique LR-pair universe
df_lr["lr_pair"] = list(
    zip(
        df_lr["ligand"],
        df_lr["receptor"]
    )
)

all_lr_pairs = set(
    df_lr["lr_pair"].unique()
)


# -----------------------------
# define disease-risk LR pairs
# -----------------------------
# A LR pair is disease-risk if either:
# 1) its ligand is a disease ligand, OR
# 2) its receptor is a disease receptor

risk_lr_pairs = set(
    (lig, rec)
    for lig, rec in all_lr_pairs
    if (
        lig in disease_ligands
        or rec in disease_receptors
    )
)


print(
    "n total LR-pair universe:",
    len(all_lr_pairs)
)

print(
    "n disease-risk LR pairs:",
    len(risk_lr_pairs)
)


# save LR-pair universe with risk annotation
lrpair_universe_df = pd.DataFrame(
    list(all_lr_pairs),
    columns=[
        "ligand",
        "receptor"
    ]
)

lrpair_universe_df[
    "ligand_is_disease"
] = lrpair_universe_df[
    "ligand"
].isin(
    disease_ligands
)

lrpair_universe_df[
    "receptor_is_disease"
] = lrpair_universe_df[
    "receptor"
].isin(
    disease_receptors
)

lrpair_universe_df[
    "is_risk_lr_pair"
] = (
    lrpair_universe_df[
        "ligand_is_disease"
    ]
    |
    lrpair_universe_df[
        "receptor_is_disease"
    ]
)


lrpair_universe_df.to_csv(
    os.path.join(
        out_path,
        f"{disease_name}_LRpair_universe_"
        f"with_either_end_risk_annotation_"
        f"{bandwidth}_mindonors{MIN_DONORS}.csv"
    ),
    index=False
)


# -----------------------------
# helper function: LR-pair enrichment
# -----------------------------
def lr_pair_enrichment_test(
    lr_pairs_in_cellpair,
    all_lr_pairs,
    risk_lr_pairs
):
    """
    Fisher enrichment test for LR pairs.

    lr_pairs_in_cellpair:
        set of unique (ligand, receptor) pairs
        in one source-target cell pair

    all_lr_pairs:
        universe of all unique LR pairs
        after donor filter

    risk_lr_pairs:
        LR pairs where either ligand or receptor
        is disease-associated
    """

    lr_pairs_in_cellpair = set(
        lr_pairs_in_cellpair
    )

    not_in_cellpair = (
        all_lr_pairs -
        lr_pairs_in_cellpair
    )

    a = len(
        lr_pairs_in_cellpair &
        risk_lr_pairs
    )

    b = len(
        lr_pairs_in_cellpair -
        risk_lr_pairs
    )

    c = len(
        not_in_cellpair &
        risk_lr_pairs
    )

    d = len(
        not_in_cellpair -
        risk_lr_pairs
    )

    n_lr_pairs_in_cellpair = a + b
    n_lr_pairs_outside_cellpair = c + d

    if (
        n_lr_pairs_in_cellpair == 0
        or n_lr_pairs_outside_cellpair == 0
    ):
        return pd.Series(
            {
                "n_lr_pairs_in_cellpair":
                    n_lr_pairs_in_cellpair,

                "n_risk_lr_pairs_in_cellpair":
                    a,

                "risk_lr_fraction":
                    np.nan,

                "n_risk_lr_pairs_outside_cellpair":
                    c,

                "n_nonrisk_lr_pairs_in_cellpair":
                    b,

                "n_nonrisk_lr_pairs_outside_cellpair":
                    d,

                "odds_ratio":
                    np.nan,

                "pval":
                    1.0
            }
        )

    odds_ratio, pval = fisher_exact(
        [
            [a, b],
            [c, d]
        ],
        alternative="greater"
    )

    return pd.Series(
        {
            "n_lr_pairs_in_cellpair":
                n_lr_pairs_in_cellpair,

            "n_risk_lr_pairs_in_cellpair":
                a,

            "risk_lr_fraction":
                a / n_lr_pairs_in_cellpair,

            "n_risk_lr_pairs_outside_cellpair":
                c,

            "n_nonrisk_lr_pairs_in_cellpair":
                b,

            "n_nonrisk_lr_pairs_outside_cellpair":
                d,

            "odds_ratio":
                odds_ratio,

            "pval":
                pval
        }
    )


# -----------------------------
# prepare LR-pair sets for each source-target cell pair
# -----------------------------
pair_keys = [
    "source",
    "target"
]

cellpair_lr_sets = (
    df_lr
    .groupby(
        pair_keys
    )["lr_pair"]
    .apply(
        lambda x: set(x)
    )
    .reset_index(
        name="lr_pair_set"
    )
)


print(
    "n source-target cell pairs tested:",
    cellpair_lr_sets.shape[0]
)


# -----------------------------
# LR-pair enrichment by source-target pair
# -----------------------------
lrpair_results = pd.concat(
    [
        cellpair_lr_sets[
            pair_keys
        ],

        cellpair_lr_sets[
            "lr_pair_set"
        ].apply(
            lambda s:
                lr_pair_enrichment_test(
                    s,
                    all_lr_pairs,
                    risk_lr_pairs
                )
        )
    ],
    axis=1
)


lrpair_results["fdr"] = multipletests(
    lrpair_results["pval"],
    method="fdr_bh"
)[1]


lrpair_results[
    "neglog10_fdr"
] = -np.log10(
    lrpair_results[
        "fdr"
    ].clip(
        lower=1e-300
    )
)


lrpair_results[
    "neglog10_pval"
] = -np.log10(
    lrpair_results[
        "pval"
    ].clip(
        lower=1e-300
    )
)


# -----------------------------
# heatmap annotation:
# significance stars only
# -----------------------------
def fdr_to_stars(fdr):

    if pd.isna(fdr):
        return ""

    elif fdr < 0.001:
        return "***"

    elif fdr < 0.01:
        return "**"

    elif fdr < 0.05:
        return "*"

    else:
        return ""


lrpair_results[
    "sig_label"
] = lrpair_results[
    "fdr"
].apply(
    fdr_to_stars
)


# Keep this for output table only,
# not for heatmap text
lrpair_results[
    "risk_count_label"
] = (
    lrpair_results[
        "n_risk_lr_pairs_in_cellpair"
    ]
    .astype(int)
    .astype(str)
    +
    "/"
    +
    lrpair_results[
        "n_lr_pairs_in_cellpair"
    ]
    .astype(int)
    .astype(str)
)


# order columns
lrpair_results = lrpair_results[
    [
        "source",
        "target",
        "n_lr_pairs_in_cellpair",
        "n_risk_lr_pairs_in_cellpair",
        "risk_lr_fraction",
        "n_nonrisk_lr_pairs_in_cellpair",
        "n_risk_lr_pairs_outside_cellpair",
        "n_nonrisk_lr_pairs_outside_cellpair",
        "odds_ratio",
        "pval",
        "fdr",
        "neglog10_pval",
        "neglog10_fdr",
        "risk_count_label",
        "sig_label"
    ]
]


# -----------------------------
# save result table
# -----------------------------
lrpair_results.to_csv(
    os.path.join(
        out_path,
        f"{disease_name}_LRpair_either_end_"
        f"enrichment_by_cellpair_"
        f"{bandwidth}_mindonors{MIN_DONORS}.csv"
    ),
    index=False
)


# also save a sorted version
lrpair_results.sort_values(
    [
        "fdr",
        "pval"
    ],
    ascending=True
).to_csv(
    os.path.join(
        out_path,
        f"{disease_name}_LRpair_either_end_"
        f"enrichment_by_cellpair_"
        f"{bandwidth}_mindonors{MIN_DONORS}_sorted.csv"
    ),
    index=False
)


print(
    "Top enriched source-target pairs:"
)

print(
    lrpair_results
    .sort_values(
        [
            "fdr",
            "pval"
        ]
    )
    .head(20)
    [
        [
            "source",
            "target",
            "n_lr_pairs_in_cellpair",
            "n_risk_lr_pairs_in_cellpair",
            "risk_lr_fraction",
            "odds_ratio",
            "pval",
            "fdr",
            "sig_label"
        ]
    ]
)


# -----------------------------
# plotting helper
# -----------------------------
def plot_lrpair_heatmap(
    res_df,
    value_col,
    annot_col,
    title,
    out_file,
    cmap="YlOrRd",
    cbar_label=""
):
    heat_data = res_df.pivot(
        index="target",
        columns="source",
        values=value_col
    )

    annot_data = res_df.pivot(
        index="target",
        columns="source",
        values=annot_col
    )

    # Keep axes in the latest
    # hd_cell_type_map.csv order
    source_order = [
        x
        for x in CELL_TYPE_ORDER
        if x in heat_data.columns
    ]

    target_order = [
        x
        for x in CELL_TYPE_ORDER
        if x in heat_data.index
    ]

    heat_data = heat_data.reindex(
        index=target_order,
        columns=source_order
    )

    annot_data = (
        annot_data
        .reindex(
            index=target_order,
            columns=source_order
        )
        .fillna("")
    )

    # Larger figure to accommodate larger fonts
    plt.figure(
        figsize=(14, 11)
    )

    ax = sns.heatmap(
        heat_data,
        annot=annot_data,
        fmt="",
        cmap=cmap,
        linewidths=0.5,
        linecolor="white",

        cbar_kws={
            "label": cbar_label
        },

        annot_kws={
            "fontsize": 15,
            "fontweight": "bold",
            "color": "black"
        }
    )

    # Force all annotation texts,
    # including *, **, ***, to be black
    for text in ax.texts:
        text.set_color("black")
        text.set_fontsize(15)
        text.set_fontweight("bold")

    # -----------------------------
    # Axis labels
    # -----------------------------
    ax.set_xlabel(
        "Source",
        fontsize=18
    )

    ax.set_ylabel(
        "Target",
        fontsize=18
    )

    # -----------------------------
    # Title
    # -----------------------------
    ax.set_title(
        title,
        fontsize=20,
        pad=18
    )

    # -----------------------------
    # Cell type labels
    # -----------------------------
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        fontsize=14
    )

    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=14
    )

    # -----------------------------
    # Colorbar
    # -----------------------------
    cbar = ax.collections[0].colorbar

    cbar.ax.tick_params(
        labelsize=14
    )

    cbar.set_label(
        cbar_label,
        fontsize=17
    )

    plt.tight_layout()

    plt.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()


# -----------------------------
# heatmap 1:
# LR-pair enrichment
# fill = -log10(FDR)
# text = FDR significance stars only
# -----------------------------
plot_lrpair_heatmap(
    lrpair_results,

    value_col="neglog10_fdr",

    annot_col="sig_label",

    title=(
        f"LR-pair Enrichment for "
        f"{disease_name} Disease-risk Interactions\n"
        f"Risk LR pair = ligand OR receptor "
        f"is disease-associated; "
        f"Donor filter: n_donors >= {MIN_DONORS}"
    ),

    out_file=os.path.join(
        plot_dir,
        f"source_target_LRpair_either_end_"
        f"enrichment_heatmap_"
        f"{disease_name}_{bandwidth}_"
        f"mindonors{MIN_DONORS}.pdf"
    ),

    cmap="YlOrRd",

    cbar_label="-log10(FDR)"
)


# -----------------------------
# heatmap 2:
# odds ratio
# fill = odds ratio
# text = FDR significance stars only
# -----------------------------
lrpair_results[
    "odds_ratio_plot"
] = lrpair_results[
    "odds_ratio"
].replace(
    [
        np.inf,
        -np.inf
    ],
    np.nan
)


plot_lrpair_heatmap(
    lrpair_results,

    value_col="odds_ratio_plot",

    annot_col="sig_label",

    title=(
        f"LR-pair Disease-risk Enrichment "
        f"Odds Ratio for {disease_name}\n"
        f"Risk LR pair = ligand OR receptor "
        f"is disease-associated; "
        f"Donor filter: n_donors >= {MIN_DONORS}"
    ),

    out_file=os.path.join(
        plot_dir,
        f"source_target_LRpair_either_end_"
        f"enrichment_OR_heatmap_"
        f"{disease_name}_{bandwidth}_"
        f"mindonors{MIN_DONORS}.png"
    ),

    cmap="YlOrRd",

    cbar_label="Odds ratio"
)


print("Done.")