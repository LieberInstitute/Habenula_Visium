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

out_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'liana'
)
os.makedirs(os.path.join(out_path, "figure", "habenula_extracellular"), exist_ok=True)

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info_split.csv')
sample_info = pd.read_csv(sample_info_path)
sample_id = sample_info['sample_id'].iloc[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]
tissue_id = sample_info['tissue_id'].iloc[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

scale_json = here(
    'processed-data', '01_spaceranger', "five_samples_10_2025", sample_id,
    'outs', 'binned_outputs', 'square_008um', 'spatial',
    'scalefactors_json.json'
)

file_path = os.path.join(out_path, f"adata/adata_extracellular_withcelltype.h5ad")
adata = sc.read(file_path)
adata = adata[adata.obs['region'] == 'habenula', :].copy()
adata.var.index.name = None

###### run separately for each sample ######
adata.obs = adata.obs.rename(columns={"tissue_section": "tissue_id"})
adata = adata[adata.obs["tissue_id"].isin([tissue_id])].copy()
adata = adata[adata.obs["region"] == "habenula"].copy()

adata = adata[:, ~adata.var_names.str.startswith("DEPRECATED_")].copy()
adata.var_names = adata.var_names.astype(str)
adata.var_names_make_unique()

###### run liana+ ######

# adata.layers['counts'] = adata.X.copy()
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# figure
sc.pl.spatial(adata, color=['cell_type'], spot_size=80, palette="tab20")
plt.savefig(os.path.join(out_path, "figure","habenula_extracellular", f"{tissue_id}_Spot_clusters.png"), dpi=300, bbox_inches='tight')
plt.close()

# Spatial Connectivity
plot, _ = li.ut.query_bandwidth(coordinates=adata.obsm['spatial'], start=0, end=200, interval_n=2)
plot
plot.save(
    os.path.join(out_path, "figure","habenula_extracellular", f"{tissue_id}_Spatial_connectivity_bandwidth.png"),
    dpi=300,
    width=6,  
    height=4,
    verbose=False
)

with open(scale_json, 'r') as f:
    scale_factors = json.load(f)

#   Number of pixels that reaches 1 Visium spot in distance, which is
#   recommended for Visium data. Here we generalize for Visium HD
good_bandwidth = 100 / scale_factors['microns_per_pixel']

# check
res = li.resource.select_resource('consensus')  # ['ligand', 'receptor']
x_name, y_name = 'ligand', 'receptor'

x_vals = set(res[x_name].astype(str).unique())
y_vals = set(res[y_name].astype(str).unique())
genes  = set(map(str, adata.var_names))

hit_lig = len(genes & x_vals)
hit_rec = len(genes & y_vals)
print(f"[PRECHECK] ligand hitting: {hit_lig}, receptor hitting: {hit_rec}")

# check numbers in liana 
print("ligand number:", len(x_vals))
print("receptor number:", len(y_vals))


li.ut.spatial_neighbors(adata, bandwidth=100, cutoff=0.1, kernel='gaussian', set_diag=True)
fig = li.pl.connectivity(adata, idx=0, size=1, figure_size=(6, 5))
fig.save(
    os.path.join(out_path, "figure", "habenula_extracellular", f"{tissue_id}_Spatial_connectivity.png"),
    dpi=300,
    width=6,   
    height=5,
    verbose=False
)


lrdata = li.mt.bivariate(adata,
                resource_name='consensus', # NOTE: uses HUMAN gene symbols!
                local_name='cosine', # Name of the function
                global_name="morans", # Name global function
                n_perms=100, # Number of permutations to calculate a p-value
                mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                add_categories=True, # Whether to add local categories to the results
                nz_prop=0.01, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True
                )

out_file = os.path.join(out_path, f"extracellular_lrdata_{tissue_id}.h5ad")
lrdata.write(out_file)

#Global Summaries
lrdata.var.sort_values("mean", ascending=False).head(5)
lrdata.var.sort_values("std", ascending=False).head(5)
lrdata.var.sort_values("morans", ascending=False).head(5)

a = lrdata.var.sort_values("morans", ascending=False).index[0]

b = lrdata.var.sort_values("std", ascending=False).index[0]

c = lrdata.var.sort_values("mean", ascending=False).index[0]

# local
# NOTE: reset params as plotnine seems to change them
sc.set_figure_params(dpi=80, dpi_save=300, format='png', frameon=False, transparent=True, figsize=[5,5])
sc.pl.spatial(
    lrdata,
    color=[a, b, c],
    spot_size=80,      
    vmax=1,
    cmap='magma',
    show=False    
)

plt.savefig(os.path.join(out_path, "figure","habenula_extracellular",f"{tissue_id}_top_Ligand-Receptor_local.png"),
            dpi=300, bbox_inches='tight')
plt.close()

# Permutation-based
sc.pl.spatial(lrdata, layer='pvals', color=[a, b, c], spot_size=80, cmap="magma_r")
plt.savefig(os.path.join(out_path, "figure","habenula_extracellular", f"{tissue_id}_top_Ligand-Receptor_permutation.png"),
            dpi=300, bbox_inches='tight')
plt.close()

# Local Categories
sc.pl.spatial(
    lrdata,
    layer='cats',
    color=[a, b, c],
    spot_size=80,      # ← manually set spot size
    cmap="coolwarm",
    show=False
)
plt.savefig(os.path.join(out_path, "figure", "habenula_extracellular", f"{tissue_id}_top_Ligand-Receptor_local_category.png"),
            dpi=300, bbox_inches='tight')
plt.close()


li.multi.nmf(lrdata, n_components=None, inplace=True, random_state=0, max_iter=200, verbose=True)

# Extract the variable loadings
lr_loadings = li.ut.get_variable_loadings(lrdata, varm_key='NMF_H').set_index('index')

# Extract the factor scores
factor_scores = li.ut.get_factor_scores(lrdata, obsm_key='NMF_W')

nmf = sc.AnnData(X=lrdata.obsm['NMF_W'],
                 obs=lrdata.obs,
                 var=pd.DataFrame(index=lr_loadings.columns),
                 uns=lrdata.uns,
                 obsm=lrdata.obsm)

sc.pl.spatial(nmf, color=[*nmf.var.index, None], spot_size=80, size=1, ncols=2, show=False)

plt.savefig(os.path.join(out_path, "figure", "habenula_extracellular", f"{tissue_id}_Intercellular_Patterns.png"),
            dpi=300, bbox_inches='tight')
plt.close()

