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

num_factors = 12
adata_in_path = here(
    'processed-data', '09_HD_cell_level', f'{sample_id}_pre_bin2cell.h5ad'
)
adata_final_path = here(
    'processed-data', '09_HD_cell_level', f'{sample_id}.h5ad'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'bin2cell_out', 
    f'{sample_id}.csv'
)
ficture_input_path = here(
    'processed-data', '10_HD_bin_level', 'ficture_harmony', 'ficture_outputs',
    'normalized', 'analysis', 'nF12.d_12', 'transcripts_joined.tsv.gz'
)
factor_cols = [
    'factor_K1', 'factor_K2', 'factor_K3', 'factor_P1', 'factor_P2', 'factor_P3'
]

os.makedirs(os.path.dirname(out_path), exist_ok=True)

#   Read in FICTURE clusters and subset to this sample
ficture_input = pd.read_csv(
    ficture_input_path,
    sep = '\t',
    usecols = ['barcode', 'sample_id'] + factor_cols
)
ficture_input = ficture_input[ficture_input['sample_id'] == sample_id]

#   Read in AnnData for this sample, but replace counts assay with zeros.
#   Instead of genes, create [num_factors] columns in this assay that will
#   correspond to "scores" for each FICTURE cluster
adata = sc.read(adata_in_path)
adata_ficture = sc.AnnData(
    X = np.zeros((adata.shape[0], num_factors), dtype = np.float32),
    obs = adata.obs,
    obsm = adata.obsm,
    uns = adata.uns
)

#   Join FICTURE cluster info into the AnnData
adata_ficture.obs[factor_cols] = (
    ficture_input
        .dropna()
        .groupby('barcode')
        .first()
        [factor_cols]
)

#   Loop through and add up scores for each potential FICTURE cluster for each
#   bin. Scores are equal to the probabilities of the top three clusters
for factor_num in range(num_factors):
    for factor_rank in range(3):
        #   Grab non-NA rows matching this particular cluster
        mask = ~adata_ficture.obs[f'factor_K{factor_rank + 1}'].isna() & \
            (adata_ficture.obs[f'factor_K{factor_rank + 1}'] == factor_num)
        
        #   Add the probability for this cluster to the corresponding column
        adata_ficture.X[mask, factor_num] += adata_ficture.obs[f'factor_P{factor_rank + 1}'][mask]

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

#   Sanity check: the number and identity of cells should be identical to the
#   ordinary bin2cell run (it should be the same input data and the algorithm
#   should be deterministic)
adata_final = sc.read(adata_final_path)
assert all(adata_final.obs.index == adata_ficture.obs.index)

#   Export a CSV containing cell key and scores for each cluster
cluster_df = pd.DataFrame(
    adata_ficture.X.toarray()
)
cluster_df.columns = [f'FICTURE_{i}' for i in range(num_factors)]
cluster_df['key'] = [f'{i}_{sample_id}' for i in adata_ficture.obs.index]
cluster_df.to_csv(out_path, index = False)

session_info.show()
