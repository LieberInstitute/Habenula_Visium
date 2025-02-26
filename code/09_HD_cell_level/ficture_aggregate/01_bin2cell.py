import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import datetime
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

# sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
# with open(sample_id_path, 'r') as f:
#     all_samples = f.read().splitlines()
# sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]
sample_id = 'H1-MVPY9BW_A1_8433'

adata_in_path = here(
    'processed-data', '09_HD_cell_level', f'{sample_id}_pre_bin2cell.h5ad'
)
adata_out_path = here(
    'processed-data', '09_HD_cell_level', 'ficture_aggregate', 'bin2cell_out', 
    f'{sample_id}.h5ad'
)
ficture_input_path = here(
    'processed-data', '09_HD_cell_level', 'ficture_aggregate', 'synthetic_input.tsv.gz'
)
mpp = 0.3

os.makedirs(os.path.dirname(adata_out_path), exist_ok=True)

adata = sc.read(adata_in_path)

ficture_input = pd.read_csv(ficture_input_path, sep = '\t')

for i in range(3):
    ficture_input[f'factor_K{i+1}'] = np.random.randint(0, 12, ficture_input.shape[0])

adata_ficture = sc.AnnData(
    X = np.zeros((adata.shape[0], 12), dtype = np.float32),
    obs = adata.obs,
    obsm = adata.obsm,
    uns = adata.uns
)
adata_ficture.obs[['factor_K1', 'factor_K2', 'factor_K3']] = (
    ficture_input
        .groupby('barcode')
        .first()
        [['factor_K1', 'factor_K2', 'factor_K3']]
)

for factor_num in range(12):
    for factor_rank in range(3):
        adata_ficture.X[
            ~adata_ficture.obs[f'factor_K{factor_rank + 1}'].isna() & 
            (adata_ficture.obs[f'factor_K{factor_rank + 1}'] == factor_num),
            factor_num
        ] += factor_rank + 1

adata_ficture.X = csr_matrix(adata_ficture.X)
adata_ficture = b2c.bin_to_cell(
    adata_ficture, labels_key="labels_joint",
    spatial_keys=["spatial", "spatial_cropped_150_buffer"]
)

sc.write(adata_out_path, adata_ficture)
session_info.show()
