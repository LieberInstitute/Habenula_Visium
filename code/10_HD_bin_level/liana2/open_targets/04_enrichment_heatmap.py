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
parser.add_argument("--bandwidth", type=float, required=True)
parser.add_argument("--disease", type=str, required=True)
args = parser.parse_args()

bandwidth = args.bandwidth
disease_name = args.disease

# For 8 donors, keep LR interactions significant in at least 4 donors
MIN_DONORS = 4

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
    "open_targets",
)

out_path = here(
    "processed-data",
    "10_HD_bin_level",
    "no_secondary",
    "liana2",
    "open_targets",
)

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(out_path, exist_ok=True)


# -----------------------------
# load disease-risk genes from OpenTargets
# -----------------------------
disease_pairs = pd.read_csv(
    os.path.join(
        out_path,
        f"{disease_name}_LorR_interactions.csv",
    )
)

# expected columns:
# genesymbol_intercell_source, genesymbol_intercell_target
disease_ligands = set(
    disease_pairs["genesymbol_intercell_source"].dropna().unique()
)

disease_receptors = set(
    disease_pairs["genesymbol_intercell_target"].dropna().unique()
)

print("n disease ligands:", len(disease_ligands))
print("n disease receptors:", len(disease_receptors))


# -----------------------------
# load LIANA interactions
# -----------------------------
df = pd.read_csv(
    f"/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/"
    f"processed-data/10_HD_bin_level/no_secondary/liana2/table/"
    f"significant_interactions_across_donors_{bandwidth}.csv"
)

# keep only needed columns
df = df[
    [
        "donor_id",
        "source",
        "target",
        "ligand",
        "receptor",
    ]
].dropna().copy()

print("input interaction rows:", df.shape[0])
print("n donors:", df["donor_id"].nunique())


# -----------------------------
# filter LR interactions by donor support
# -----------------------------
donor_support = (
    df.groupby(
        [
            "source",
            "target",
            "ligand",
            "receptor",
        ]
    )["donor_id"]
    .nunique()
    .reset_index(name="n_donors")
)

print(
    "unique source-target-LR before donor filter:",
    donor_support.shape[0],
)

df = df.merge(
    donor_support,
    on=[
        "source",
        "target",
        "ligand",
        "receptor",
    ],
    how="left",
)

df = df[df["n_donors"] >= MIN_DONORS].copy()

print(
    f"interaction rows after donor filter n_donors >= {MIN_DONORS}:",
    df.shape[0],
)

print(
    f"unique source-target-LR after donor filter n_donors >= {MIN_DONORS}:",
    df[
        [
            "source",
            "target",
            "ligand",
            "receptor",
        ]
    ]
    .drop_duplicates()
    .shape[0],
)


# -----------------------------
# define universes
# -----------------------------
all_ligands = set(df["ligand"].unique())
all_receptors = set(df["receptor"].unique())

# restrict disease genes to genes actually observable in the interaction universe
disease_ligands = disease_ligands & all_ligands
disease_receptors = disease_receptors & all_receptors

print("n disease ligands in universe:", len(disease_ligands))
print("n disease receptors in universe:", len(disease_receptors))
print("n total ligand universe:", len(all_ligands))
print("n total receptor universe:", len(all_receptors))


# -----------------------------
# helper function
# -----------------------------
def enrichment_test(genes_in_pair, all_genes, disease_genes):
    """
    genes_in_pair:
        set of unique genes in a source-target pair

    all_genes:
        universe set

    disease_genes:
        disease-associated gene set within universe
    """

    not_in_pair = all_genes - genes_in_pair

    a = len(genes_in_pair & disease_genes)
    b = len(genes_in_pair - disease_genes)
    c = len(not_in_pair & disease_genes)
    d = len(not_in_pair - disease_genes)

    # edge case
    if (a + b == 0) or (c + d == 0):
        return pd.Series(
            {
                "n_genes_in_pair": len(genes_in_pair),
                "n_disease_in_pair": a,
                "odds_ratio": float("nan"),
                "pval": 1.0,
            }
        )

    odds_ratio, pval = fisher_exact(
        [
            [a, b],
            [c, d],
        ],
        alternative="greater",
    )

    return pd.Series(
        {
            "n_genes_in_pair": len(genes_in_pair),
            "n_disease_in_pair": a,
            "odds_ratio": odds_ratio,
            "pval": pval,
        }
    )


# -----------------------------
# prepare unique gene sets for each source-target pair
# -----------------------------
pair_keys = ["source", "target"]

ligand_sets = (
    df.groupby(pair_keys)["ligand"]
    .apply(lambda x: set(x.unique()))
    .reset_index(name="gene_set")
)

receptor_sets = (
    df.groupby(pair_keys)["receptor"]
    .apply(lambda x: set(x.unique()))
    .reset_index(name="gene_set")
)


# -----------------------------
# ligand enrichment
# -----------------------------
ligand_results = ligand_sets.copy()

ligand_results = pd.concat(
    [
        ligand_results[pair_keys],
        ligand_results["gene_set"].apply(
            lambda s: enrichment_test(
                s,
                all_ligands,
                disease_ligands,
            )
        ),
    ],
    axis=1,
)

ligand_results["fdr"] = multipletests(
    ligand_results["pval"],
    method="fdr_bh",
)[1]

ligand_results["neglog10_fdr"] = -np.log10(
    ligand_results["fdr"].clip(lower=1e-300)
)


# -----------------------------
# receptor enrichment
# -----------------------------
receptor_results = receptor_sets.copy()

receptor_results = pd.concat(
    [
        receptor_results[pair_keys],
        receptor_results["gene_set"].apply(
            lambda s: enrichment_test(
                s,
                all_receptors,
                disease_receptors,
            )
        ),
    ],
    axis=1,
)

receptor_results["fdr"] = multipletests(
    receptor_results["pval"],
    method="fdr_bh",
)[1]

receptor_results["neglog10_fdr"] = -np.log10(
    receptor_results["fdr"].clip(lower=1e-300)
)


# -----------------------------
# save result tables
# -----------------------------
ligand_results.to_csv(
    os.path.join(
        out_path,
        f"{disease_name}_ligand_enrichment_by_cellpair_{bandwidth}.csv",
    ),
    index=False,
)

receptor_results.to_csv(
    os.path.join(
        out_path,
        f"{disease_name}_receptor_enrichment_by_cellpair_{bandwidth}.csv",
    ),
    index=False,
)


# -----------------------------
# plotting helper
# -----------------------------
def plot_heatmap(
    res_df,
    value_col,
    annot_col,
    title,
    out_file,
    cmap="YlOrRd",
):
    heat_data = res_df.pivot(
        index="target",
        columns="source",
        values=value_col,
    )

    annot_data = res_df.pivot(
        index="target",
        columns="source",
        values=annot_col,
    )

    # Slightly larger figure for larger fonts
    plt.figure(figsize=(12, 10))

    ax = sns.heatmap(
        heat_data,
        annot=annot_data.round(3),
        fmt="",
        cmap=cmap,
        linewidths=0.5,
        linecolor="white",
        annot_kws={
            "fontsize": 12,
        },
        cbar_kws={
            "label": "-log10(FDR)",
        },
    )

    # Axis labels
    ax.set_xlabel(
        "Source",
        fontsize=16,
    )

    ax.set_ylabel(
        "Target",
        fontsize=16,
    )

    # Title
    ax.set_title(
        title,
        fontsize=18,
        pad=15,
    )

    # X-axis tick labels
    ax.set_xticklabels(
        ax.get_xticklabels(),
        rotation=45,
        ha="right",
        fontsize=13,
    )

    # Y-axis tick labels
    ax.set_yticklabels(
        ax.get_yticklabels(),
        rotation=0,
        fontsize=13,
    )

    # Colorbar formatting
    cbar = ax.collections[0].colorbar

    cbar.ax.tick_params(
        labelsize=13,
    )

    cbar.set_label(
        "-log10(FDR)",
        fontsize=15,
    )

    plt.tight_layout()

    plt.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close()


# -----------------------------
# heatmap 1: ligand enrichment
# fill = -log10(FDR), annotation = FDR
# -----------------------------
plot_heatmap(
    ligand_results,
    value_col="neglog10_fdr",
    annot_col="fdr",
    title=(
        f"Ligand Enrichment for {disease_name} Disease Genes\n"
        f"Fill = -log10(FDR), Text = FDR"
    ),
    out_file=os.path.join(
        plot_dir,
        f"source_target_ligand_enrichment_heatmap_{disease_name}_{bandwidth}.png",
    ),
    cmap="YlOrRd",
)


# -----------------------------
# heatmap 2: receptor enrichment
# fill = -log10(FDR), annotation = FDR
# -----------------------------
plot_heatmap(
    receptor_results,
    value_col="neglog10_fdr",
    annot_col="fdr",
    title=(
        f"Receptor Enrichment for {disease_name} Disease Genes\n"
        f"Fill = -log10(FDR), Text = FDR"
    ),
    out_file=os.path.join(
        plot_dir,
        f"source_target_receptor_enrichment_heatmap_{disease_name}_{bandwidth}.png",
    ),
    cmap="YlGnBu",
)

print("Done.")