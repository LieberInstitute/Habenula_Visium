import scanpy as sc
import decoupler as dc
import plotnine as p9
import liana as li
import numpy as np
import pandas as pd
import os
from pyhere import here
import matplotlib.pyplot as plt
import anndata as ad

out_path = here(
    'processed-data', '10_HD_bin_level', 'LIANA'
)

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()
    sample_id = all_samples[1 - 1]
    sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

file_path = os.path.join(out_path, f"lrdata_{sample_id}.h5ad")
lrdata = ad.read_h5ad(file_path)

###### run LIANA+ ######
# figure
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

