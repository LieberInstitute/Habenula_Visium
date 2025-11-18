import scanpy as sc
import decoupler as dc
import plotnine as p9
import liana as li
import numpy as np
import pandas as pd
import os
from pyhere import here
import matplotlib.pyplot as plt

out_path = here(
    'processed-data', '10_HD_bin_level', 'LIANA'
)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()
    sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

file_path = os.path.join(out_path, f"adata/adata_{sample_id}.h5ad")
adata = sc.read(file_path)

###### run LIANA+ ######

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# figure
sc.pl.spatial(adata, color=[None, 'cell_type'], size=30, palette="tab20")
plt.savefig(os.path.join(out_path, "figure", f"{sample_id}_Spot_clusters.png"), dpi=300, bbox_inches='tight')
plt.close()

# Spatial Connectivity
plot, _ = li.ut.query_bandwidth(coordinates=adata.obsm['spatial'], start=0, end=500, interval_n=20)
plot
plot.save(
    os.path.join(out_path, "figure", f"{sample_id}_Spatial_connectivity_bandwidth.png"),
    dpi=300,
    width=6,  # 单位 inch，可自调
    height=4,
    verbose=False
)

li.ut.spatial_neighbors(adata, bandwidth=80, cutoff=0.1, kernel='gaussian', set_diag=True)
fig = li.pl.connectivity(adata, idx=0, size=1.3, figure_size=(6, 5))
fig.save(
    os.path.join(out_path, "figure", f"{sample_id}_Spatial_connectivity.png"),
    dpi=300,
    width=6,   # 单位 inch，可按 figure_size 填
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
                nz_prop=0.001, # Minimum expr. proportion for ligands/receptors and their subunits
                use_raw=False,
                verbose=True
                )

out_file = os.path.join(out_path, f"lrdata_{sample_id}.h5ad")
lrdata.write(out_file)

#Global Summaries
lrdata.var.sort_values("mean", ascending=False).head(5)
lrdata.var.sort_values("std", ascending=False).head(5)
lrdata.var.sort_values("morans", ascending=False).head(5)

# local
# NOTE: reset params as plotnine seems to change them
sc.set_figure_params(dpi=80, dpi_save=300, format='png', frameon=False, transparent=True, figsize=[5,5])
sc.pl.spatial(
    lrdata,
    color=['APP^GPC1', 'PSAP^CELSR1'],
    size=10,
    vmax=1,
    alpha_img=0.7,
    cmap='magma',
    show=False    # 不弹出窗口
)
plt.savefig(os.path.join(out_path, "figure", f"{sample_id}_top_Ligand-Receptor_local.png"),
            dpi=300, bbox_inches='tight')
plt.close()

# Permutation-based
sc.pl.spatial(lrdata, layer='pvals', color=['APP^GPC1', 'PSAP^CELSR1'], size=10, alpha_img=0.7, cmap="magma_r")
plt.savefig(os.path.join(out_path, "figure", f"{sample_id}_top_Ligand-Receptor_permutation.png"),
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

sc.pl.spatial(nmf, color=[*nmf.var.index, None], size=10, alpha_img=0.7, ncols=2)
plt.savefig(os.path.join(out_path, "figure", f"{sample_id}_Intercellular_Patterns.png"),
            dpi=300, bbox_inches='tight')
plt.close()

lr_loadings.sort_values("Factor2", ascending=False).head(10)

