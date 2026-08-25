"""
Prepare gene sets from trios.csv for enrichment analysis.

Reads trios.csv, filters to 'Intersect' trios, and creates a long-format
gene_sets CSV with TFs and genes as separate gene sets per cell type.

Output: trios_intersect_gene_sets.csv (columns: gene_set_name, gene)
"""

import pandas as pd
import os
from pyhere import here

# -----------------------------
# paths
# -----------------------------
ot_dir = str(
    here(
        "code", "10_HD_bin_level", "liana2", "open_targets"
    )
)

trios = pd.read_csv(os.path.join(ot_dir, "trios.csv"))

# Filter to level 2 trios
trios["test_level"] = trios["test_level"].astype(str)
trios = trios[trios["test_level"] == "2"].copy()
print(f"n Intersect trios: {len(trios)}")
print(f"cell types: {sorted(trios['cell_type'].unique())}")

# Build gene sets: TFs and genes per cell type
records = []

for ct in trios["cell_type"].unique():
    ct_df = trios[trios["cell_type"] == ct]

    # TF gene set
    tfs = ct_df["TF"].dropna().unique()
    for tf in tfs:
        records.append({"gene_set_name": f"{ct}_TFs", "gene": tf})

    # Gene gene set
    genes = ct_df["gene"].dropna().unique()
    for gene in genes:
        records.append({"gene_set_name": f"{ct}_genes", "gene": gene})

gene_sets_df = pd.DataFrame(records).drop_duplicates()

out_file = os.path.join(ot_dir, "trios_intersect_gene_sets.csv")
gene_sets_df.to_csv(out_file, index=False)

print(f"\nOutput: {out_file}")
print(f"n gene sets: {gene_sets_df['gene_set_name'].nunique()}")
print(gene_sets_df.groupby("gene_set_name").size().to_string())
