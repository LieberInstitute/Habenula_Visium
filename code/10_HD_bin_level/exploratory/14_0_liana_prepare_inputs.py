import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import datetime
import anndata as ad

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()[:3]

ad_in_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    '{}.h5ad'
)
out_dir = here('processed-data', '10_HD_bin_level', 'LIANA', 'adata')

adata_list = []
for sample_id in all_samples:
    adata = sc.read(str(ad_in_paths).format(sample_id))
    adata.obs['sample_id'] = sample_id
    adata.obs['key'] = adata.obs.index + '_' + adata.obs['sample_id']
    adata.obs.index = adata.obs['key']
    adata_list.append(adata)

adata = ad.concat(adata_list, axis=0)

#   Perform filtering and log normalization similar to that done for the
#   cell-level object in R
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#   Filter out the artifact in H1-MVPY9BW_A1_8433
adata = adata[
    (adata.obs['sample_id'] != 'H1-MVPY9BW_A1_8433') |
    (adata.obsm['spatial'][:,0] <= 34223),
    :
]

sc.write(os.path.join(out_dir, 'cellular.h5ad'), adata)
