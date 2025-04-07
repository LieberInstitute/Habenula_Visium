import os
from pyhere import here
import session_info
import bin2cell as b2c
import datetime
import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import sys
import gzip
from scipy.io import mmread
import HERGAST
import torch

sample_id = 'H1-W369TJK_D1_9090'
sr_dir = here(
    'processed-data', '01_spaceranger', sample_id, 'outs', 'binned_outputs',
    'square_008um'
)
sr_spatial_dir = here(
    'processed-data', '01_spaceranger', sample_id, 'outs', 'spatial'
)
out_path = here('processed-data', '10_HD_bin_level', 'hergast_adata.h5ad')
plot_dir = here('plots', '10_HD_bin_level')
raw_image_path = here('raw-data', 'images', 'vis-hd', f'{sample_id}.tif')
hvg_path = here('processed-data', '09_HD_cell_level', 'HVGs.txt')
random_seed = 0

os.makedirs(plot_dir, exist_ok=True)
plt.rcParams["figure.figsize"] = (15, 10)

#   Report if a GPU is available
if torch.cuda.is_available():
    print("Running on GPU.")
else:
    print("Running on CPU only!")

################################################################################
#   Build and preprocess AnnData
################################################################################

print(f"{datetime.datetime.now()} | Building and preprocessing AnnData")

#   Read in spaceranger outputs into an AnnData
adata = b2c.read_visium(
    sr_dir,
    source_image_path = raw_image_path,
    spaceranger_image_path = sr_spatial_dir
)

#   Use Ensembl IDs for var_names
adata.var_names = adata.var['gene_ids']
adata.var_names.name = None

#   Take in-tissue bins containing some expressed genes, and genes expressed in
#   at least one bin
adata = adata[adata.obs['in_tissue'] == 1, :]
adata = adata[
    np.sum(adata.X, axis = 1) > 0, np.sum(adata.X, axis = 0) > 0
]

#   Read in HVGs, which we'll subset to during PCA
with open(hvg_path, 'r') as f:
    hvg = f.read().split('\n')[: -1]

adata.var['is_hvg'] = False
adata.var.loc[hvg, 'is_hvg'] = True

#   Normalize expression and compute 200 PCs
sc.pp.normalize_total(adata, target_sum = 1, exclude_highly_expressed = True)
sc.pp.scale(adata)
sc.pp.pca(adata, n_comps = 200, random_state = random_seed, mask_var = 'is_hvg')

################################################################################
#  Run HERGAST
################################################################################

#   Construct relational graph
print(f"{datetime.datetime.now()} | Constructing relational graph")
HERGAST.utils.Cal_Spatial_Net(adata, k_cutoff = 8)
HERGAST.utils.Cal_Expression_Net(adata, dim_reduce = 'PCA')

#   The authors recommend dividing the data into batches with 10,000 to 20,000
#   spots each. Calculate the number of segments in each dimension to achieve
#   this
num_splits = int(np.round((adata.shape[0] / 10000) ** 0.5))
assert num_splits > 1

#   Train model
print(f"{datetime.datetime.now()} | Training model using num_batch_x_y=({num_splits}, {num_splits})")
train_HERGAST = HERGAST.Train_HERGAST(
    adata, batch_data = True, num_batch_x_y = (num_splits, num_splits),
    spatial_net_arg = {'k_cutoff': 8, 'verbose': False},
    exp_net_arg = {'verbose': False}, dim_reduction = 'PCA', device_idx = 0
)
train_HERGAST.train_HERGAST(n_epochs = 200)

#   Perform final clustering
print(f"{datetime.datetime.now()} | Performing final clustering")
sc.pp.neighbors(adata, use_rep='HERGAST')
sc.tl.umap(adata)
sc.tl.leiden(adata, random_state = random_seed, resolution = 0.3)

sc.pl.spatial(
    adata, color = 'leiden', title = 'HERGAST Clusters', palette = 'tab20',
    show = False
)
plt.savefig(os.path.join(plot_dir, f'{sample_id}_HERGAST.png'))
plt.close('all')

sc.write(out_path, adata)

session_info.show()
