# Cell type–specific analysis: Group ligand–receptor pairs by cell type and select top pairs for each.
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
from scipy import sparse
from pathlib import Path

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))

#   Read input files
in_dir = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana'
)
plot_dir= here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana', 'figure','habenula'
)
os.makedirs(plot_dir, exist_ok=True)

in_files = [ os.path.join(in_dir, f) for f in os.listdir(in_dir) if re.compile(r'.*\.h5ad$').match(f) ]

base_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure/habenula"
if task_id == 1:
    in_files = [f for f in in_files if "extracellular" not in f and "lrdata" in f]
    output_dir1 = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_files"
    output_dir2 = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_specific_interactions"
    plot_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure/habenula/celltype_specific_interactions"
    global_dir = os.path.join(base_dir, "spatial_top_pairs_global")
    specific_dir = os.path.join(base_dir, "spatial_top_pairs_specific")
else:
    in_files = [f for f in in_files if "extracellular" in f and "lrdata" in f]
    output_dir1 = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_files_extracellular"
    output_dir2 = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/celltype_specific_interactions_extracellular"
    plot_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure/habenula_extracellular/celltype_specific_interactions_extracellular"
    base_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure/habenula_extracellular"
    global_dir = os.path.join(base_dir, "spatial_top_pairs_global_extracellular")
    specific_dir = os.path.join(base_dir, "spatial_top_pairs_specific_extracellular")

os.makedirs(output_dir1, exist_ok=True)
os.makedirs(output_dir2, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(base_dir, exist_ok=True)
os.makedirs(global_dir, exist_ok=True)
os.makedirs(specific_dir, exist_ok=True)
# ===================================================================

# -------- helper: pick an (n_obs x n_vars) matrix from AnnData --------
def pick_matrix(adata, prefer=None):
    """
    Return a 2D (n_obs x n_vars) matrix for interactions and its source name.
    Priority: user 'prefer' -> .X (if width==n_vars) -> 'cats' -> 'pvals'.
    """
    if prefer is not None:
        if prefer in adata.layers and adata.layers[prefer].shape[1] == adata.n_vars:
            return adata.layers[prefer], prefer
        raise ValueError(f"Requested layer '{prefer}' not found or wrong shape.")
    if adata.X is not None and getattr(adata.X, "shape", (0, 0))[1] == adata.n_vars:
        return adata.X, "X"
    for name in ("cats", "pvals"):
        if name in adata.layers and adata.layers[name].shape[1] == adata.n_vars:
            return adata.layers[name], name
    raise ValueError(
        f"No (n_obs x n_vars) matrix found. .X is {getattr(adata.X,'shape',None)}; "
        f"available layers: {list(adata.layers.keys())}"
    )

# -------- helper: compute per-cell_type mean scores for one AnnData --------
def celltype_means_for_donor(adata, donor_id):
    # choose matrix
    M, source = pick_matrix(adata)
    # do NOT densify if big: compute group means sparsely
    ct = adata.obs['cell_type'].values
    var_names = adata.var_names
    # unique cell types in this donor
    ct_unique = pd.Index(pd.Series(ct, dtype='category').cat.categories)
    if len(ct_unique) == 0:
        ct_unique = pd.Index(np.unique(ct))
    rows = []
    # compute mean per cell_type without densifying whole matrix
    for c in ct_unique:
        idx = np.where(ct == c)[0]
        if len(idx) == 0:
            continue
        sub = M[idx]  # submatrix (n_idx x n_vars)
        # mean across rows for this cell type
        if sparse.issparse(sub):
            m = np.array(sub.mean(axis=0)).ravel()
        else:
            m = sub.mean(axis=0)
            if hasattr(m, "A"):  # numpy.matrix
                m = m.A.ravel()
        s = pd.Series(m, index=var_names, name=c)
        rows.append(s)
    wide = pd.DataFrame(rows)
    wide.index.name = 'cell_type'
    wide.reset_index(inplace=True)
    # long format
    long = wide.melt(id_vars='cell_type', var_name='interaction', value_name='mean_score')
    long['donor_id'] = donor_id
    return long, wide, source

files = [in_files[0], in_files[1], in_files[2]]
donor_results = []
sources_used = []
var_name_sets = []
for f in files:
    donor_id = Path(f).stem.replace('lrdata_', '').replace('extracellular_', '')
    adata = sc.read(f)
    # compute
    long, wide, src = celltype_means_for_donor(adata, donor_id)
    donor_results.append((donor_id, long, wide))
    sources_used.append(src)
    var_name_sets.append(set(adata.var_names))

print("Matrix sources used (in order):", sources_used)
#Matrix sources used (in order): ['X', 'X', 'X']

# -------------------------------------------------------
# get ranking of ligand-receptor pairs by cell type and donor
donor_sorted = {}

for donor_id, long_df, _ in donor_results:
    sorted_df = (
        long_df
        .sort_values(['cell_type', 'mean_score'], ascending=[True, False])
        .groupby('cell_type', group_keys=False)
        .apply(lambda x: x.reset_index(drop=True))
    )
    donor_sorted[donor_id] = sorted_df

# --------------------------------------------------------
# rank ligand-receptor pairs by mean score across donors

# combine all long dataframes across donors
combined_long = pd.concat([x[1] for x in donor_results], ignore_index=True)

# calculate average mean_score for each (cell_type, interaction) pair across donors
avg_df = (
    combined_long
    .groupby(['cell_type', 'interaction'], as_index=False)
    ['mean_score']
    .mean()
)

# rank within each cell_type
avg_sorted = (
    avg_df
    .sort_values(['cell_type', 'mean_score'], ascending=[True, False])
    .groupby('cell_type', group_keys=False)
    .apply(lambda x: x.reset_index(drop=True))
)

# select top 10 pairs per cell type by average mean_score
avg_top10 = (
    avg_df
    .sort_values(['cell_type', 'mean_score'], ascending=[True, False])
    .groupby('cell_type', group_keys=False)
    .head(10)
    .reset_index(drop=True)
)
avg_top10.to_csv(
    os.path.join(output_dir1, "celltype_top10_by_mean_score.csv"),
    index=False
)

# get the union of top10 interactions across cell types
top10_interactions = avg_top10["interaction"].drop_duplicates()
avg_top10_all_ct = avg_df[avg_df["interaction"].isin(top10_interactions)].copy()

heatmap_mat = avg_top10_all_ct.pivot_table(
    index="interaction",
    columns="cell_type",
    values="mean_score",
    aggfunc="max"
)

heatmap_mat.to_csv(
    os.path.join(output_dir1, "celltype_top10_union_all_celltypes_matrix.csv")
)

# =======================================================
# for each cell_type, create a file with interactions ranked by mean_score per donor
cell_types = combined_long["cell_type"].unique()

# output directory
os.makedirs(output_dir1, exist_ok=True)

# loop over cell types
for ct in cell_types:
    # get subset for this cell type
    subset = combined_long.query("cell_type == @ct")
    # pivot：raw=interaction，col=donor_id，value=mean_score
    pivot_df = subset.pivot_table(
        index="interaction", columns="donor_id", values="mean_score"
    )
    pivot_df = pivot_df.dropna(how="any")  # drop rows with any NaN
    # add mean across donors
    pivot_df["Mean"] = pivot_df.mean(axis=1)    
    # rank by mean
    pivot_df = pivot_df.sort_values("Mean", ascending=False)
    safe_ct = re.sub(r'[\\/:"*?<>|]+', "_", str(ct))
    # output file path
    out_path = os.path.join(output_dir1, f"{safe_ct}_interactions.csv")
    pivot_df.to_csv(out_path)

# =======================================================

# 1) calculate mean score per (cell_type, interaction) across donors
ct_int_means = (
    combined_long
    .groupby(['cell_type', 'interaction'], as_index=False)['mean_score']
    .mean()
    .rename(columns={'mean_score': 'ct_mean'})
)

# 2) summarize within cell types: n, mean, std
g = ct_int_means.groupby('interaction')
tmp = g['ct_mean'].agg(n_ct='size', g_mean='mean', g_std='std').reset_index()
df = ct_int_means.merge(tmp, on='interaction', how='left')

# 3) compute leave-one-out mean and std for each cell type
eps = 1e-12
df['others_mean'] = np.where(
    df['n_ct'] > 1,
    (df['g_mean'] * df['n_ct'] - df['ct_mean']) / (df['n_ct'] - 1),
    np.nan
)

S_total   = (df['g_std']**2) * (df['n_ct'] - 1)              # Σ(x - mean)^2
d         = df['ct_mean'] - df['g_mean']
S_others  = S_total - (df['n_ct'] / (df['n_ct'] - 1).replace(0, np.nan)) * (d**2)
denom     = (df['n_ct'] - 2)
others_var = np.where(denom > 0, S_others / denom, np.nan)
others_var = np.clip(others_var, 0, None)

df['others_std'] = np.sqrt(others_var + eps)

# 4) define specificity metrics relative to "other cell types"
df['z_spec'] = (df['ct_mean'] - df['others_mean']) / (df['others_std'] + eps)
df['abs_z']  = df['z_spec'].abs()

df_sorted = (
    df.sort_values(['cell_type', 'abs_z'], ascending=[True, False])
      .reset_index(drop=True)
)

# 5) select top pairs per cell type by absolute z_spec
top10_per_ct = (
    df.sort_values(['cell_type', 'abs_z'], ascending=[True, False])
      .groupby('cell_type', group_keys=False)
      .head(10)
      .reset_index(drop=True)
)[['cell_type','interaction','ct_mean','others_mean','others_std','z_spec','abs_z']]

# 6) save the top pairs per cell type
os.makedirs(output_dir2, exist_ok=True)

df_sorted.to_csv(os.path.join(output_dir2, "celltype_specific_interactions_all.csv"), index=False)
top10_per_ct.to_csv(os.path.join(output_dir2, "celltype_specific_interactions_top10.csv"), index=False)

# =======================================================
# Plot top pairs across donors
os.makedirs(plot_dir, exist_ok=True)

import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, Normalize
import numpy as np
import re, os

# Define color maps
cmap_div = plt.get_cmap('coolwarm')   # For mixed positive and negative values (diverging)
cmap_warm = plt.get_cmap('Reds')      # For all positive values (warm tones)
cmap_cool = plt.get_cmap('Blues')     # For all negative values (cool tones)

for ct, sub in top10_per_ct.groupby('cell_type'):
    # Sort by z_spec in descending order
    sub_sorted = sub.sort_values('z_spec', ascending=False).copy()
    zmin, zmax = sub_sorted['z_spec'].min(), sub_sorted['z_spec'].max()
    # 🔹 Automatically choose color scheme based on value distribution
    if zmin < 0 and zmax > 0:
        # Both positive and negative values → use diverging colormap (centered at 0)
        norm = TwoSlopeNorm(vmin=zmin, vcenter=0, vmax=zmax)
        cmap = cmap_div
    elif zmax <= 0:
        # All values are negative → use cool colors (darker = more negative)
        norm = Normalize(vmin=zmin, vmax=zmax)
        cmap = cmap_cool
    else:
        # All values are positive → use warm colors (darker = more positive)
        norm = Normalize(vmin=zmin, vmax=zmax)
        cmap = cmap_warm
    # Map each z_spec value to a color
    colors = cmap(norm(sub_sorted['z_spec'].values))
    # 🔸 Plot
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.barh(sub_sorted['interaction'], sub_sorted['z_spec'], color=colors)
    ax.axvline(0, color='gray', lw=1)  # Vertical line at 0 for reference
    ax.set_xlabel("Specificity Z-score (z_spec)")
    ax.set_title(f"{ct}: Top 10 Specific Interactions (descending)")
    ax.invert_yaxis()  # Show highest scores at the top
    # 🔹 Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("z_spec")
    # 🔹 Save the plot safely
    safe_ct = re.sub(r'[\\/:"*?<>|]+', "_", str(ct))
    out_path = os.path.join(plot_dir, f"{safe_ct}_top10_specific.png")
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()

# =======================================================
# Visualize top ligand-receptor pairs in spatial context

os.makedirs(global_dir, exist_ok=True)
os.makedirs(specific_dir, exist_ok=True)

# ==== utils ====
def safe_name(x: str) -> str:
    return re.sub(r'[\\/:"*?<>|]+', "_", str(x))

def get_vals_for_pairs(adata, pairs):
    """return (n_spots, n_pairs), missing pairs are skipped."""
    valid = [p for p in pairs if p in adata.var_names]
    if not valid:
        return None, []
    X = adata.X
    idx = adata.var_names.get_indexer(valid)
    vals = X[:, idx]
    if hasattr(vals, "toarray"):
        vals = vals.toarray()
    return np.asarray(vals), valid

def get_global_minmax_for_pairs(h5ads, pairs):
    gmin, gmax = np.inf, -np.inf
    for f in h5ads:
        ad = sc.read(f)
        ad.obsm['spatial'] = np.asarray(ad.obsm['spatial'])
        vals, valid = get_vals_for_pairs(ad, pairs)
        if vals is None:
            continue
        gmin = min(gmin, np.nanmin(vals))
        gmax = max(gmax, np.nanmax(vals))
    if not np.isfinite(gmin) or not np.isfinite(gmax):
        gmin, gmax = 0.0, 1.0
    return float(gmin), float(gmax)

# ==== params ====
top_n = 10
# ------------------------------------------------------
# A) Global Top-10 pairs across samples & cell types
# ------------------------------------------------------
global_pair_mean = (
    combined_long
    .groupby('interaction', as_index=False)['mean_score']
    .mean()
    .rename(columns={'mean_score':'global_mean'})
)

global_rank = (
    global_pair_mean
    .sort_values('global_mean', ascending=False)
)

global_rank.to_csv("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/global_interaction_ranking_mean.csv")

global_top = (
    global_pair_mean
    .sort_values('global_mean', ascending=False)
    .head(top_n)['interaction'].tolist()
)

print("Global Top pairs:", global_top)

# all samples together to get unified color scale
gmin, gmax = get_global_minmax_for_pairs(in_files, global_top)

for f in in_files:
    ad = sc.read(f)
    sample_id = os.path.basename(f).replace(".h5ad", "").replace("lrdata_", "")
    ad.obsm['spatial'] = np.asarray(ad.obsm['spatial'])
    # columns existing in this sample
    valid = [p for p in global_top if p in ad.var_names]
    if not valid:
        print(f"[skip] {sample_id}: none of global_top in var_names.")
        continue    # sele
    fig = sc.pl.spatial(
        ad,
        color=valid,
        cmap='viridis',
        vmin=gmin, vmax=gmax,
        ncols=2,
        spot_size=80,
        show=False,
        return_fig=True
    )
    out = os.path.join(global_dir, f"global_top{top_n}_{safe_name(sample_id)}.pdf")
    fig.savefig(out, bbox_inches='tight')
    plt.close(fig)
    print("[saved]", out)

# ------------------------------------------------------
# B) Cell-type specific Top-10 pairs (from top10_per_ct)
# ------------------------------------------------------
# Assume you already have `top10_per_ct` (each row = cell_type × interaction × metrics)
for ct, sub in top10_per_ct.groupby('cell_type'):
    # Select top_n interactions for this cell type
    ct_pairs = (
        sub.sort_values('abs_z', ascending=False)
           .head(top_n)['interaction'].tolist()
    )
    if not ct_pairs:
        continue
    # Compute a unified color scale across samples for this set of pairs
    cmin, cmax = get_global_minmax_for_pairs(in_files, ct_pairs)    
    # Loop through each donor/sample
    for f in in_files:
        ad = sc.read(f)
        sample_id = os.path.basename(f).replace(".h5ad", "").replace("lrdata_", "")
        ad.obsm['spatial'] = np.asarray(ad.obsm['spatial'])
        # Check which pairs are present in this sample
        valid_pairs = [p for p in ct_pairs if p in ad.var_names]
        if not valid_pairs:
            print(f"[skip] {ct} | {sample_id}: none of ct_pairs found.")
            continue
        fig = sc.pl.spatial(
            ad,
            color=valid_pairs,
            cmap='viridis',
            vmin=cmin, vmax=cmax,
            ncols=2,
            spot_size=80,
            show=False,
            return_fig=True
        )
        # Save output figure
        out = os.path.join(specific_dir, f"{safe_name(ct)}_top{top_n}_{safe_name(sample_id)}.pdf")
        fig.savefig(out, bbox_inches='tight')
        plt.close(fig)
        print("[saved]", out)

session_info.show()
