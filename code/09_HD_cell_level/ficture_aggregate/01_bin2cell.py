import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

#   Grab sample ID using the array task ID
sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()
sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

adata_in_path = here(
    'processed-data', '09_HD_cell_level', f'{sample_id}_pre_bin2cell.h5ad'
)
adata_final_path = here('processed-data', '09_HD_cell_level', f'{sample_id}.h5ad')
out_path = here(
    'processed-data', '09_HD_cell_level', 'ficture_aggregate', 'bin2cell_out', 
    f'{sample_id}.csv'
)
ficture_input_path = here(
    'processed-data', '09_HD_cell_level', 'ficture_aggregate', 'synthetic_input.tsv.gz'
)

os.makedirs(os.path.dirname(out_path), exist_ok=True)

ficture_input = pd.read_csv(ficture_input_path, sep = '\t')

#   For now, randomly assign clusters to each bin just to have some working data
for i in range(3):
    ficture_input[f'factor_K{i+1}'] = np.random.randint(0, 12, ficture_input.shape[0])

#   Read in AnnData for this sample, but replace counts assay with zeros.
#   Instead of genes, create 12 columns in this assay that will correspond
#   to "scores" for each FICTURE cluster
adata = sc.read(adata_in_path)
adata_ficture = sc.AnnData(
    X = np.zeros((adata.shape[0], 12), dtype = np.float32),
    obs = adata.obs,
    obsm = adata.obsm,
    uns = adata.uns
)

#   Join FICTURE cluster info into the AnnData
adata_ficture.obs[['factor_K1', 'factor_K2', 'factor_K3']] = (
    ficture_input
        .groupby('barcode')
        .first()
        [['factor_K1', 'factor_K2', 'factor_K3']]
)

#   Loop through and add up scores for each potential FICTURE cluster for each
#   bin. Add 3 points for the top-ranked cluster call (i.e. 'factor_K1'), 2
#   for 'factor_K2', and 1 for 'factor_K3'. In this way, each bin that was
#   assigned a FICTURE cluster (some will be null) has a total of 6 points to
#   allocate across 3 clusters
for factor_num in range(12):
    for factor_rank in range(3):
        adata_ficture.X[
            ~adata_ficture.obs[f'factor_K{factor_rank + 1}'].isna() & 
            (adata_ficture.obs[f'factor_K{factor_rank + 1}'] == factor_num),
            factor_num
        ] += 3 - factor_rank

#   Apply bin_to_cell on the scores matrix. The idea is that a cell will add
#   up scores for each potential cluster across constituent bins. A cell by this
#   definition is a sort of distribution across potentially multiple FICTURE
#   clusters, but of course it can be placed back into a discrete category by
#   assigning the top-scoring cluster later
adata_ficture.X = csr_matrix(adata_ficture.X)
adata_ficture = b2c.bin_to_cell(
    adata_ficture, labels_key="labels_joint",
    spatial_keys=["spatial", "spatial_cropped_150_buffer"]
)

#   Sanity check: the number of cells should be identical to the ordinary
#   bin2cell run (it should be the same input data and the algorithm should be
#   deterministic)
adata_final = sc.read(adata_final_path)
assert adata_final.shape[0] == adata_ficture.shape[0], f'Actual vs. expected shapes: {adata_ficture.shape[0]} vs. {adata_final.shape[0]}'

#   Export a CSV containing cell key and scores for each cluster
cluster_df = pd.DataFrame(
    adata_ficture.X.toarray()
)
cluster_df.columns = [f'FICTURE_{i}' for i in range(12)]
cluster_df['key'] = [f'{i}_{sample_id}' for i in cluster_df.index]
cluster_df.to_csv(out_path, index = False)

session_info.show()
