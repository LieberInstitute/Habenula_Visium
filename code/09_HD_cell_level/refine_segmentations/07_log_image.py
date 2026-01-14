
#   Create log-scaled gene-expression images. In many cases in the original
#   images, the poor dynamic range results in many primary segmentations not
#   showing obvious gene expression

import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import pandas as pd

task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = pd.read_csv(sample_info_path)
sample_id = sample_info.iloc[task_id]['sample_id']

stardist_dir = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'stardist'
)
pre_out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    f'{sample_id}_pre_bin2cell.h5ad'
)
mpp = 0.3

################################################################################
#   Write the image
################################################################################

adata = sc.read_h5ad(pre_out_path)
adata.obs['sum_umi'] = adata.X.sum(axis = 1)

#   Create an image from gene counts
b2c.grid_image(
    adata,
    "sum_umi",
    mpp=mpp,
    log1p = True,
    sigma=5,
    save_path=os.path.join(
        stardist_dir, f'gex_log_{sample_id}.tiff'
    )
)

session_info.show()
