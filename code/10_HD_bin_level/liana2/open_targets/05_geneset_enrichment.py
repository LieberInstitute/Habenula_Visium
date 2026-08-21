"""
Gene set enrichment against OpenTargets disease-risk genes.

For each gene set of interest, tests whether it is enriched for disease-risk
genes using Fisher's exact test (one-sided, greater). Applies BH FDR correction
across all tests.

Inputs:
    --gene_sets: CSV with columns 'gene_set_name' and 'gene' (long format)
    --disease: Disease name(s) matching {disease}_risk_genes_01thr.csv files
               produced by 01_prep_open_targets_tsv.R. Can specify multiple.
    --background: Single-column CSV (or one gene per line) of background genes.
    --output_dir: Directory for output files.

Outputs:
    - CSV with enrichment results (one row per gene_set × disease combination)
    - Summary dot plot (gene set × disease) if multiple diseases provided
"""

import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from pyhere import here


# -----------------------------
# args
# -----------------------------
parser = argparse.ArgumentParser(
    description="Test gene sets for enrichment in OpenTargets disease-risk genes."
)
parser.add_argument(
    "--gene_sets",
    required=True,
    help="CSV with columns 'gene_set_name' and 'gene'.",
)
parser.add_argument(
    "--disease",
    required=True,
    nargs="+",
    help="Disease name(s) matching {disease}_risk_genes_01thr.csv files.",
)
parser.add_argument(
    "--background",
    required=True,
    help="Single-column CSV or text file listing background genes.",
)
parser.add_argument(
    "--output_dir",
    default=None,
    help="Output directory. Defaults to processed-data/.../open_targets/.",
)
parser.add_argument(
    "--risk_gene_dir",
    default=None,
    help=(
        "Directory containing {disease}_risk_genes_01thr.csv files. "
        "Defaults to processed-data/10_HD_bin_level/no_secondary/liana2/open_targets/."
    ),
)

args = parser.parse_args()

# -----------------------------
# paths
# -----------------------------
default_ot_dir = str(
    here(
        "processed-data", "10_HD_bin_level", "no_secondary", "liana2", "open_targets"
    )
)

risk_gene_dir = args.risk_gene_dir if args.risk_gene_dir else default_ot_dir
output_dir = args.output_dir if args.output_dir else default_ot_dir
os.makedirs(output_dir, exist_ok=True)

# -----------------------------
# load data
# -----------------------------
# Gene sets of interest
gene_sets_df = pd.read_csv(args.gene_sets)
assert {"gene_set_name", "gene"}.issubset(
    gene_sets_df.columns
), "gene_sets file must have columns 'gene_set_name' and 'gene'"
gene_sets_df = gene_sets_df[["gene_set_name", "gene"]].dropna().drop_duplicates()

# Background genes
bg_df = pd.read_csv(args.background, header=0)
# Handle single-column file regardless of column name
background = set(bg_df.iloc[:, 0].dropna().astype(str).unique())

print(f"n background genes: {len(background)}")
print(f"n gene sets: {gene_sets_df['gene_set_name'].nunique()}")

# -----------------------------
# load disease-risk genes for each disease
# -----------------------------
disease_risk = {}
for disease in args.disease:
    path = os.path.join(risk_gene_dir, f"{disease}_risk_genes_01thr.csv")
    df_risk = pd.read_csv(path)
    # Expect column 'gene'
    if "gene" not in df_risk.columns:
        # Try first column
        df_risk = df_risk.rename(columns={df_risk.columns[0]: "gene"})
    risk_genes = set(df_risk["gene"].dropna().astype(str).unique())
    # Restrict to background
    risk_genes = risk_genes & background
    disease_risk[disease] = risk_genes
    print(f"disease={disease}: {len(risk_genes)} risk genes in background")

# -----------------------------
# restrict gene sets to background
# -----------------------------
gene_sets_df = gene_sets_df[gene_sets_df["gene"].isin(background)].copy()

print(
    f"gene sets after restricting to background: "
    f"{gene_sets_df['gene_set_name'].nunique()} sets, "
    f"{len(gene_sets_df)} gene-set entries"
)

# -----------------------------
# enrichment tests
# -----------------------------
results = []

for disease, risk_genes in disease_risk.items():
    for gs_name, group in gene_sets_df.groupby("gene_set_name"):
        gs_genes = set(group["gene"].unique())

        a = len(gs_genes & risk_genes)
        b = len(gs_genes - risk_genes)
        c = len(risk_genes - gs_genes)
        d = len(background - gs_genes - risk_genes)

        if (a + b) == 0:
            odds_ratio, pval = np.nan, 1.0
        else:
            odds_ratio, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

        overlap_genes = sorted(gs_genes & risk_genes)

        results.append(
            {
                "gene_set_name": gs_name,
                "disease": disease,
                "n_genes_in_set": a + b,
                "n_risk_genes": len(risk_genes),
                "n_overlap": a,
                "odds_ratio": odds_ratio,
                "pval": pval,
                "overlapping_genes": ";".join(overlap_genes),
            }
        )

results_df = pd.DataFrame(results)

# FDR correction across all tests jointly
if len(results_df) > 0:
    results_df["fdr"] = multipletests(results_df["pval"], method="fdr_bh")[1]
else:
    results_df["fdr"] = []

results_df["neglog10_fdr"] = -np.log10(results_df["fdr"].clip(lower=1e-300))
results_df["neglog10_pval"] = -np.log10(results_df["pval"].clip(lower=1e-300))

# Sort
results_df = results_df.sort_values(["fdr", "pval"]).reset_index(drop=True)

# -----------------------------
# save
# -----------------------------
diseases_label = "_".join(args.disease)
out_file = os.path.join(
    output_dir, f"geneset_enrichment_{diseases_label}.csv"
)
results_df.to_csv(out_file, index=False)
print(f"Results saved to: {out_file}")
print(f"n significant (FDR < 0.05): {(results_df['fdr'] < 0.05).sum()}")
print(results_df.head(20).to_string(index=False))

print("Done.")
