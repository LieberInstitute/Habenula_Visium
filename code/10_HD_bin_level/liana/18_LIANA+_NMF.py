# NMF
# Global NMF: Combine the three samples and perform global NMF analysis.
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
import anndata as ad

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))

#   Read input files
in_dir = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana'
)

in_files = [ os.path.join(in_dir, f) for f in os.listdir(in_dir) if re.compile(r'.*\.h5ad$').match(f) ]

if task_id == 1:
    in_files = [f for f in in_files if "extracellular" not in f]
    plot_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure/habenula/NMF"
    table_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF"
    data_suffix = ""
    out_file = os.path.join(in_dir, f"merged_five_files.h5ad")
else:
    in_files = [f for f in in_files if "extracellular" in f and "lrdata" in f]
    plot_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/figure/habenula_extracellular/NMF"
    table_dir = "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/liana/table/NMF_extracellular"
    data_suffix = ""
    out_file = os.path.join(in_dir, f"merged_five_files_extracellular.h5ad")

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(table_dir, exist_ok=True)

# read the three files
adatas = [sc.read_h5ad(f) for f in in_files]

# check the dimensions of each object
for i, a in enumerate(adatas):
    print(f"{os.path.basename(in_files[i])}: {a.shape}")

# merge the five objects
adata_merged = ad.concat(adatas, join="inner", label="batch", keys=[os.path.basename(f) for f in in_files])

# save the merged object
# adata_merged.write(out_file)

print(f"Merged AnnData saved to {out_file}")

# Shift spatial coordinates to avoid overlap
# get unique batches
batches = adata_merged.obs["batch"].unique()

x_offset = 0
gap = 300  # gap between samples

for b in batches:
    idx = adata_merged.obs["batch"] == b
    coords = adata_merged.obsm["spatial"][idx, :].copy()
    coords[:, 0] = coords[:, 0] - coords[:, 0].min()
    width = coords[:, 0].max()
    coords[:, 0] = coords[:, 0] + x_offset
    adata_merged.obsm["spatial"][idx, :] = coords
    x_offset += width + gap

print("✅ Spatial coordinates shifted successfully.")

# Perform NMF analysis on the merged data
#   Read in DataFrames of ligand-receptor stats for each donor and concatenate
li.multi.nmf(adata_merged, n_components=None, inplace=True, random_state=0, max_iter=300, verbose=True)

# Extract the variable loadings
lr_loadings = li.ut.get_variable_loadings(adata_merged, varm_key='NMF_H').set_index('index')

lr_loadings.to_csv(os.path.join(table_dir, f"NMF_H_loadings{data_suffix}.csv"))

# Extract the factor scores
factor_scores = li.ut.get_factor_scores(adata_merged, obsm_key='NMF_W')

nmf = sc.AnnData(X=adata_merged.obsm['NMF_W'],
                 obs=adata_merged.obs,
                 var=pd.DataFrame(index=lr_loadings.columns),
                 uns=adata_merged.uns,
                 obsm=adata_merged.obsm)

sc.pl.spatial(nmf, color=[*nmf.var.index, None], spot_size=80, size=1, ncols=2, show=False)

plt.savefig(os.path.join(plot_dir, f"overall_Intercellular_Patterns{data_suffix}.png"),
            dpi=300, bbox_inches='tight')
plt.close()

# =====================================
# generate a pair list
pairs = adata_merged.var.index.tolist()
ligands = [p.split("^")[0] for p in pairs]
receptors = [p.split("^")[1] for p in pairs]
all_genes = sorted(set(ligands + receptors))
pd.Series(all_genes, name="gene").to_csv(os.path.join(table_dir, f"universe_genes2{data_suffix}.txt"), index=False, header=False)

adata_merged_outer = ad.concat(adatas, join="outer", label="batch", keys=[os.path.basename(f) for f in in_files])
pairs = adata_merged_outer.var.index.tolist()
ligands = [p.split("^")[0] for p in pairs]
receptors = [p.split("^")[1] for p in pairs]
all_genes = sorted(set(ligands + receptors))
pd.Series(all_genes, name="gene").to_csv(os.path.join(table_dir, f"universe_genes{data_suffix}.txt"), index=False, header=False)

# ===================================== 
# Heatmap: Create a heatmap with factors on the y-axis, cell types on the x-axis, and fill color as the average factor scores.
# Annotate factor scores with cell types
cell_type_df = adata_merged.obs[['cell_type']].reset_index()

# Merge factor scores with cell type annotations
factor_scores_annot = factor_scores.merge(cell_type_df, on='key', how='left')

# Calculate average factor scores for each cell type
factor_cols = [c for c in factor_scores_annot.columns if c.startswith("Factor")]

factor_means = (
    factor_scores_annot
    .groupby("cell_type")[factor_cols]
    .mean()
    .reset_index()
)

# make heatmap
heatmap_data = factor_means.set_index("cell_type")

import seaborn as sns

plt.figure(figsize=(6, 6))  
sns.heatmap(
    heatmap_data,
    cmap="Oranges",        
    annot=True,            
    fmt=".2f",
    cbar_kws={"label": "Average Factor Score"}
)

plt.xlabel("Factor")
plt.ylabel("Cell Type")
plt.title("Average Factor Scores by Cell Type")
plt.yticks(rotation=0)     
plt.tight_layout()

plt.savefig(os.path.join(plot_dir, f"celltype_factor_heatmap{data_suffix}.pdf"), dpi=300, bbox_inches="tight")

# =============================================================
# Cell type specific NMF: Perform NMF analysis for each cell type separately.
cell_types = adata_merged.obs['cell_type'].unique()

for cell_type in cell_types:
    adata_ct = adata_merged[adata_merged.obs['cell_type'] == cell_type].copy()
    if adata_ct.n_obs < 10:
        print(f"Skipping cell type {cell_type} due to insufficient observations ({adata_ct.n_obs})")
        continue
    print(f"Performing NMF for cell type: {cell_type} with {adata_ct.n_obs} observations")
    li.multi.nmf(adata_ct, n_components=None, inplace=True, random_state=0, max_iter=300, verbose=True)    
    # Extract the variable loadings
    lr_loadings_ct = li.ut.get_variable_loadings(adata_ct, varm_key='NMF_H').set_index('index')
    safe_ct = re.sub(r'[\\/:"*?<>|]+', "_", str(cell_type)) 
    lr_loadings_ct.to_csv(os.path.join(table_dir, f"NMF_H_loadings_{safe_ct}{data_suffix}.csv"))
    
    # Extract the factor scores
    factor_scores_ct = li.ut.get_factor_scores(adata_ct, obsm_key='NMF_W')
    
    nmf_ct = sc.AnnData(X=adata_ct.obsm['NMF_W'],
                        obs=adata_ct.obs,
                        var=pd.DataFrame(index=lr_loadings_ct.columns),
                        uns=adata_ct.uns,
                        obsm=adata_ct.obsm)
    
    sc.pl.spatial(nmf_ct, color=[*nmf_ct.var.index, None], spot_size=80, size=1, ncols=2, show=False)
    
    plt.savefig(os.path.join(plot_dir, f"{safe_ct}_Intercellular_Patterns{data_suffix}.png"),
                dpi=300, bbox_inches='tight')
    plt.close()

# for cell_type in ["OPC"]:
#     adata_ct = adata_merged[adata_merged.obs['cell_type'] == cell_type].copy()
#     if adata_ct.n_obs < 10:
#         print(f"Skipping cell type {cell_type} due to insufficient observations ({adata_ct.n_obs})")
#         continue
#     print(f"Performing NMF for cell type: {cell_type} with {adata_ct.n_obs} observations")
#     li.multi.nmf(adata_ct, n_components=5, inplace=True, random_state=0, max_iter=200, verbose=True)    
#     # Extract the variable loadings
#     lr_loadings_ct = li.ut.get_variable_loadings(adata_ct, varm_key='NMF_H').set_index('index')
#     safe_ct = re.sub(r'[\\/:"*?<>|]+', "_", str(cell_type)) 
#     lr_loadings_ct.to_csv(os.path.join(table_dir, f"NMF_H_loadings_{safe_ct}{data_suffix}.csv"))
    
#     # Extract the factor scores
#     factor_scores_ct = li.ut.get_factor_scores(adata_ct, obsm_key='NMF_W')
    
#     nmf_ct = sc.AnnData(X=adata_ct.obsm['NMF_W'],
#                         obs=adata_ct.obs,
#                         var=pd.DataFrame(index=lr_loadings_ct.columns),
#                         uns=adata_ct.uns,
#                         obsm=adata_ct.obsm)
    
#     sc.pl.spatial(nmf_ct, color=[*nmf_ct.var.index, None], spot_size=80, size=1, ncols=2, show=False)
    
#     plt.savefig(os.path.join(plot_dir, f"{safe_ct}_Intercellular_Patterns{data_suffix}.png"),
#                 dpi=300, bbox_inches='tight')
#     plt.close()


# ====================================
# Bar plots: For each factor, create a horizontal bar plot showing the top 10 ligand-receptor pairs with the highest loadings.
import math

sns.set_theme(style="whitegrid", context="talk")

factors = list(lr_loadings.columns)
n = len(factors)
ncols = 2
nrows = math.ceil(n / ncols)

fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4*nrows), squeeze=False)
plt.subplots_adjust(hspace=0.7, wspace=0.3)

palette = sns.color_palette("Blues_d", 10)

for i, factor in enumerate(factors):
    ax = axes[i // ncols][i % ncols]
    top10 = lr_loadings[factor].nlargest(10).sort_values(ascending=True)  
    ax.barh(top10.index, top10.values, color=palette)
    ax.set_title(f"{factor}: Top 10 LR pairs", fontsize=14, weight="bold")
    ax.set_xlabel("Loading")
    ax.set_ylabel("Ligand^Receptor")
    ax.tick_params(axis="y", labelsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

if nrows*ncols > n:
    for j in range(n, nrows*ncols):
        fig.delaxes(axes[j // ncols][j % ncols])

fig.tight_layout(pad=1.5)
fig.savefig(os.path.join(plot_dir, f"NMF_top10_pairs_per_factor{data_suffix}.pdf"), bbox_inches="tight", dpi=300)
plt.close(fig)
print("✅ Saved: NMF_top10_pairs_per_factor.pdf")

# ====================================
# bar plots for cell type specific NMF



